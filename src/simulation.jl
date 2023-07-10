const DUMP_FILE_HEADER = "\n------ n ------ j ------ ω ------ q --- p_wf --- p_wh " *
    "--- T_wh --- T_us --- q_ch ------ f ----- ∂f ------ λ --- q_nr --- cnvg"

function simulate_well(; kwargs...)
    well = rand_well(; kwargs...)
    data = simulate_well(well; kwargs...)
    return well, data
end

simulate_well(well; kwargs...) = simulate_well(well, rand_choke_opening(; kwargs...); kwargs...)

function simulate_well(
    well, choke_opening;
    max_wb_seg_len=MAX_WB_SEG_LEN, max_fl_seg_len=MAX_FL_SEG_LEN,
    kwargs...
)
    L_wb = well.wellbore.length # ft
    L_fl = well.flowline.length # ft

    N_t  = length(choke_opening)          # -, number of time steps
    N_wb = ceil(Int, L_wb/max_wb_seg_len) # -, number of wellbore segments
    N_fl = ceil(Int, L_fl/max_fl_seg_len) # -, number of flowline segments

    data = SimulationData(N_t, N_wb, N_fl)

    simulate_well!(
        data, well, choke_opening;
        max_wb_seg_len, max_fl_seg_len,
        kwargs...
    )

    return data
end

function simulate_well!(
    data, well, choke_opening;
    time_step=TIME_STEP, max_iter=MAX_ITER,
    rate_abs_tol=RATE_ABS_TOL, pres_abs_tol=PRES_ABS_TOL, temp_abs_tol=TEMP_ABS_TOL,
    max_wb_seg_len=MAX_WB_SEG_LEN, max_fl_seg_len=MAX_FL_SEG_LEN,
    dump_file=nothing,
    _...
)
    #* timers
    t_total = -time() # s, function timer
    t_setup = 0.0     # s, setup timer
    t_init  = 0.0     # s, initialization timer
    t_sim   = 0.0     # s, simulation timer
    t_rsrv  = 0.0     # s, reservoir calculation timer
    t_wb    = 0.0     # s, wellbore calculation timer
    t_fl    = 0.0     # s, flowline calculation timer
    t_ch    = 0.0     # s, choke calculation timer

    #* setup
    t_setup -= time() # s, setup timer

    # constant factors
    p_sc  = STD_COND_PRES     # psia, standard condition pressure
    T_sc  = STD_COND_TEMP     # °F, standard condition temperature
    ρ_wsc = STD_COND_WTR_DENS # lb/ft3, standard condition water density
    g     = GRAV_ACC          # ft/s2, gravitational acceleration
    nz    = NEAR_ZERO         # -, near zero

    # conversion factors
    α_L = LEN_CONV_FACT   # in/ft, length conversion factor
    α_t = TIME_CONV_FACT  # s/hr, time conversion factor
    α_g = GRAV_CONV_FACT  # lbm-ft/lbf-s2, gravity conversion factor
    α_q = RATE_CONV_FACT  # BPD/ft3/s, flow rate conversion factor
    α_μ = VISC_CONV_FACT  # cP/psi-s, viscosity conversion factor
    α_k = PERM_CONV_FACT  # mD/ft2, permeability conversion factor
    α_V = VALVE_CONV_FACT # BPD/gpm, valve flow rate conversion factor

    # simulation parameters
    Δt        = time_step      # hr, time step
    i_max     = max_iter       # -, maximum number of convergence iterations
    ε_q       = rate_abs_tol   # BPD, absolute tolerance on rate
    ε_p       = pres_abs_tol   # psi, absolute tolerance on pressure
    ε_T       = temp_abs_tol   # °F, absolute tolerance on temperature
    Δx_wb_max = max_wb_seg_len # ft, maximum wellbore segment length
    Δx_fl_max = max_fl_seg_len # ft, maximum flowline segment length

    # well parameters
    γ_api = well.oil.gravity         # °API
    μ_sc  = well.oil.viscosity       # cP
    c_o   = well.oil.compressibility # 1/psi
    β     = well.oil.expansion       # 1/°F
    C_p   = well.oil.capacity        # Btu/lb-°F

    p_i = well.reservoir.pressure        # psia
    k   = well.reservoir.permeability    # mD
    φ   = well.reservoir.porosity        # -
    h   = well.reservoir.thickness       # ft
    r_w = well.reservoir.radius          # ft
    c_f = well.reservoir.compressibility # 1/psi
    T_r = well.reservoir.temperature     # °F
    s   = well.reservoir.skin            # -
    
    L_wb = well.wellbore.length    # ft
    D_wb = well.wellbore.diameter  # in
    U_wb = well.wellbore.transfer  # Btu/hr-ft2-°F
    ε_wb = well.wellbore.roughness # ft

    L_fl = well.flowline.length   # ft
    D_fl = well.flowline.diameter # in
    U_fl = well.flowline.transfer # Btu/hr-ft2-°F

    T̄_a  = well.ambient.average # °F
    ΔT_a = well.ambient.annual  # °F
    ΔT_d = well.ambient.daily   # °F
    θ_T  = well.ambient.shift   # -

    p_nw = well.network.pressure # psia

    boundary = well.boundary
    choke    = well.choke

    # segment data and vectors
    p_wb_mat    = data.wellbore.pressure    # psia
    T_wb_mat    = data.wellbore.temperature # °F
    Δp_f_wb_mat = data.wellbore.friction    # psi
    T_fl_mat    = data.flowline.temperature # °F
    ψ           = data.auxilliary.transient # -

    # steps
    N_t  = length(choke_opening)     # -, number of time steps
    N_wb = ceil(Int, L_wb/Δx_wb_max) # -, number of wellbore segments
    N_fl = ceil(Int, L_fl/Δx_fl_max) # -, number of flowline segments

    # dimensions checks
    @assert length(data.time) >= N_t + 1
    @assert size(data.wellbore.pressure, 1) >= N_t + 1
    @assert size(data.wellbore.pressure, 2) >= N_wb + 1
    @assert size(data.flowline.temperature, 1) >= N_t + 1
    @assert size(data.flowline.temperature, 2) >= N_fl + 1

    # helper functions
    fvf(p, T) = exp(β*(T - T_sc) - c_o*(p - p_sc)) # RB/STB, formation volume factor

    amb_temp(n) = T̄_a + ΔT_a*cos(2*π/31536000*α_t*Δt*n + θ_T) + ΔT_d*cos(2*π/86400*α_t*Δt*n + θ_T)

    function friction_factor(N_Re)
        return if N_Re <= nz
            0.0
        elseif N_Re <= 2000
            64/N_Re
        elseif N_Re >= 4000
            β_s = α_L*ε_wb/D_wb
            (-2*log10(β_s/3.7065 - 5.0452/N_Re*log10(β_s^1.1098/2.8257 + 5.8506/N_Re^0.8981)))^-2
        else
            β_s    = α_L*ε_wb/D_wb
            f_lami = 64/2000
            f_turb = (-2*log10(β_s/3.7065 - 5.0452/4000*log10(β_s^1.1098/2.8257 + 5.8506/4000^0.8981)))^-2
            f_lami + (N_Re - 2000)*(f_turb - f_lami)/2000
        end
    end

    # constant parameters
    ρ_sc  = ρ_wsc*141.5/(γ_api + 131.5)   # lb/ft3, standard condition density
    c_t   = c_f + c_o                     # 1/psi, total compressibility
    B_r   = fvf(p_i, T_r)                 # RB/STB, reservoir FVF
    μ_r   = μ_sc/B_r                      # cP, reservoir viscosity
    H     = 4*α_t*α_μ*k/(α_k*φ*μ_r*c_t)   # ft2/hr, hydraulic diffusivity
    J     = 2*π*α_q*α_μ*k*h/(α_k*B_r*μ_r) # BPD/psi, productivity index
    g_e   = (T_r - T̄_a)/L_wb              # °F/ft, geothermal gradient
    Δx_wb = L_wb/N_wb                     # ft, wellbore segment length
    Δx_fl = L_fl/N_fl                     # ft, flowline segment length

    # dimensionless well pressure vector loop
    for n in 1:N_t
        t    = n*Δt # hr, time
        t_D  = α_t*α_μ*k*t/(α_k*φ*μ_r*c_t*r_w^2) # -, dimensionless time
        ψ[n] = dimensionless_well_pressure(boundary, t_D, s, r_w) # -, dimensionless well pressure
    end

    # additional constant parameters
    ψ_1   = ψ[1]   # -, transient function for one time step
    ∂p_wf = -ψ_1/J # psi/BPD, bottomhole pressure gradient

    # convergence parameters
    cnvg = false # convergence flag
    j    = 0     # -, number of convergence iterations
    i_sub_perf_max = div(i_max, 10) # -, maximum number of sub-performance iterations tolerated

    t_setup += time() # s, setup timer

    #* initialization
    t_init -= time() # s, initialization timer

    # temperature
    T_a = amb_temp(0)     # °F, ambient temperature
    T_fl_mat[0, :] .= T_a # °F, flowline temperature

    for i in 0:N_wb
        T_wb_mat[0, i] = T_r - g_e*i*Δx_wb # °F, wellbore temperature
    end
    
    # pressure
    p_wb_mat[0, 0] = p_i # psia, wellbore pressure

    # convergence loop
    while j < i_max && !cnvg
        j += 1 # -, number of convergence iterations

        p_wh_prior = p_wb_mat[0, N_wb] # psia, prior wellhead pressure

        # wellbore loop
        for i in 1:N_wb
            p_prev = p_wb_mat[0, i - 1] # psia, previous pressure
            T_prev = T_wb_mat[0, i - 1] # °F, previous temperature
            p      = p_wb_mat[0, i]     # psia, pressure
            T      = T_wb_mat[0, i]     # °F, temperature

            p̄ = 0.5*(p + p_prev) # psia, average pressure
            T̄ = 0.5*(T + T_prev) # °F, average temperature
            B = fvf(p̄, T̄)        # RB/STB, FVF
            ρ = ρ_sc/B           # lb/ft3, density
 
            p = p_prev - ρ*g*Δx_wb/(α_L^2*α_g) # psia, pressure

            p_wb_mat[0, i] = p # psia, wellbore pressure
        end

        p_wh = p_wb_mat[0, N_wb] # psia, wellhead pressure

        cnvg = abs(p_wh - p_wh_prior) < ε_p # convergence check
    end

    Δp_f_wb_mat[0, :] .= 0 # psi, wellbore friction pressure drop

    # initial state
    data.time[0]                      = 0                 # hr, time
    data.choke[0]                     = 0                 # -, choke opening
    data.rate[0]                      = 0                 # BPD, rate
    data.bottomhole.pressure[0]       = p_i               # psia
    data.bottomhole.temperature[0]    = T_r               # °F
    data.wellhead.pressure[0]         = p_wb_mat[0, N_wb] # psia
    data.wellhead.temperature[0]      = T̄_a               # °F
    data.upstream.pressure[0]         = p_wb_mat[0, N_wb] # psia
    data.upstream.temperature[0]      = T_a               # °F
    data.ambient.temperature[0]       = T_a               # °F
    data.auxilliary.convergence[0]    = true
    data.auxilliary.iterations[0]     = 0
    data.auxilliary.subperformance[0] = false

    t_init += time() # s, initialization timer

    #* simulation
    t_sim -= time() # s, simulation timer

    # time loop
    for n in 1:N_t
        # constant parameters
        t   = n*Δt                       # hr, time
        ω   = choke_opening[n]           # -, choke opening
        C_V = flow_coefficient(choke, ω) # gpm/sqrt(psi), valve flow coefficient
        T_a = amb_temp(n)                # °F, ambient temperature

        # sum of previous drawdowns initialization
        t_rsrv -= time() # s, reservoir calculation timer

        Δp_Σ = 0.0 # psi, sum of previous drawdowns

        # sum of previous drawdowns loop
        for m in 1:(n - 1)
            Δp_Σ += (data.rate[m] - data.rate[m - 1])*ψ[n - m + 1] # psi, sum of previous drawdowns
        end

        Δp_Σ /= J # psi, sum of previous drawdowns

        t_rsrv += time() # s, reservoir calculation timer

        # convergence loop initialization
        q_prev = data.rate[n - 1]      # BPD, previous rate
        q      = q_prev                # BPD, rate
        p_wf   = p_wb_mat[n - 1, 0]    # psia, bottomhole pressure
        p_wh   = p_wb_mat[n - 1, N_wb] # psia, wellhead pressure
        T_wh   = T_wb_mat[n - 1, N_wb] # °F, wellhead temperature
        T_us   = T_fl_mat[n - 1, N_fl] # °F, upstream temperature
        cnvg   = false                 # convergence flag
        j      = 0                     # -, number of convergence iterations

        i_sub_perf = 0     # -, number of sub-performant iterations
        sub_perf   = false # sub-performance flag

        # vector initialization
        for i in 0:N_wb
            p_wb_mat[n, i] = p_wb_mat[n - 1, i] # psia, wellbore pressure
            T_wb_mat[n, i] = T_wb_mat[n - 1, i] # °F, wellbore temperature
        end

        for i in 0:N_fl
            T_fl_mat[n, i] = T_fl_mat[n - 1, i] # °F, flowline temperature
        end

        !isnothing(dump_file) && print(dump_file, DUMP_FILE_HEADER)

        # convergence loop
        while j < i_max && !cnvg
            j += 1 # -, number of convergence iterations

            !isnothing(dump_file) && @printf dump_file "\n%8d %8d %8.6f %8.2f" n j ω q
            
            # prior values
            p_wf_prior = p_wf # psia, bottomhole pressure
            p_wh_prior = p_wh # psia, wellhead pressure
            T_wh_prior = T_wh # °F, wellhead temperature
            T_us_prior = T_us # °F, upstream temperature

            # sub-performant convergence check
            i_sub_perf += q <= nz                      # -, number of sub-performant iterations
            sub_perf    = i_sub_perf >= i_sub_perf_max # sub-performance flag

            # set production rate to zero if convergence is sub-performant
            if sub_perf
                q = 0.0
            end

            #* reservoir
            t_rsrv -= time() # s, reservoir calculation timer

            p_wf = p_i - Δp_Σ - (q - q_prev)*ψ_1/J # psia, bottomhole pressure

            t_rsrv += time() # s, reservoir calculation timer

            !isnothing(dump_file) && @printf dump_file " %8.2f" p_wf

            # ensure bottomhole pressure is greater than network pressure
            if !sub_perf && p_wf < p_nw
                !isnothing(dump_file) && @printf dump_file " BELOW NETWORK PRESSURE"
                q /= 2 # BPD, rate
                continue
            end

            #* wellbore
            t_wb -= time() # s, wellbore calculation timer

            p_wb_mat[n, 0] = p_wf  # psia, wellbore pressure
            ∂p_wh          = ∂p_wf # psi/BPD, wellhead pressure gradient
            invalid_pres   = false # invalid wellbore flag

            # wellbore loop
            for i in 1:N_wb
                p_prev = p_wb_mat[n, i - 1] # psia, previous pressure
                T_prev = T_wb_mat[n, i - 1] # psia, previous temperature
                T_old  = T_wb_mat[n - 1, i] # °F, old temperature
                p      = p_wb_mat[n, i]     # psia, pressure
                T      = T_wb_mat[n, i]     # °F, temperature
                T_e    = T_r - g_e*i*Δx_wb  # °F, geothermal temperature

                p̄ = 0.5*(p + p_prev) # psia, average pressure
                T̄ = 0.5*(T + T_prev) # °F, average temperature
                B = fvf(p̄, T̄)        # RB/STB, FVF
                ρ = ρ_sc/B           # lb/ft3, density
                μ = μ_sc/B           # cP, viscosity

                Δp_h = ρ*g*Δx_wb/(α_L^2*α_g) # psi, hydrostatic pressure drop

                u    = 4*α_L^2*q*B/(π*α_q*D_wb^2) # ft/s, velocity
                N_Re = α_μ*ρ*u*D_wb/(α_L^3*α_g*μ) # -, reynolds number
                f    = friction_factor(N_Re)      # -, friction factor

                Δp_f = 0.5*f*ρ*u^2*Δx_wb/(α_L*α_g*D_wb) # psi, friction pressure drop

                p = p_prev - Δp_h - Δp_f # psia, pressure

                invalid_pres = i < N_wb && p < p_nw
                if !sub_perf && invalid_pres
                    break
                end

                ∂u     = 4*α_L^2*B/(π*α_q*D_wb^2)      # ft/s/BPD, velocity gradient
                ∂p_wh -= f*ρ*u*∂u*Δx_wb/(α_L*α_g*D_wb) # psi/BPD, wellhead pressure gradient

                τ = 0.25*α_t*C_p*ρ*D_wb/(α_L*U_wb) # s, relaxation time
                T = (T_old/(α_t*Δt) + u*T_prev/Δx_wb + T_e/τ)/(1/(α_t*Δt) + u/Δx_wb + 1/τ) # °F

                p_wb_mat[n, i]    = p    # psia, wellbore pressure
                T_wb_mat[n, i]    = T    # °F, wellbore temperature
                Δp_f_wb_mat[n, i] = Δp_f # psi, wellbore friction pressure drop
            end

            t_wb += time() # s, wellbore calculation timer

            if !sub_perf && invalid_pres
                !isnothing(dump_file) && @printf dump_file " WELLBORE BELOW NETWORK PRESSURE"
                q /= 2 # BPD, rate
                continue
            end
            
            p_wh = p_wb_mat[n, N_wb] # psia, wellhead pressure
            T_wh = T_wb_mat[n, N_wb] # °F, wellhead temperature

            !isnothing(dump_file) && @printf dump_file " %8.2f %8.2f" p_wh T_wh

            # ensure wellbore pressure is greater than network pressure
            if !sub_perf && p_wh < p_nw
                !isnothing(dump_file) && @printf dump_file " WELLHEAD BELOW NETWORK PRESSURE"
                q /= 2 # BPD, rate
                continue
            end

            #* flowline
            t_fl -= time() # s, flowline calculation timer

            T_fl_mat[n, 0] = T_wh # °F, flowline temperature

            # flowline loop
            for i in 1:N_fl
                T_prev = T_fl_mat[n, i - 1] # °F, previous temperature
                T_old  = T_fl_mat[n - 1, i] # °F, old temperature
                T      = T_fl_mat[n, i]     # °F, temperature

                T̄ = 0.5*(T + T_prev) # °F, average temperature
                B = fvf(p_wh, T̄)     # RB/STB, FVF
                ρ = ρ_sc/B           # lb/ft3, density

                u = 4*α_L^2*q*B/(π*α_q*D_fl^2)     # ft/s, velocity
                τ = 0.25*α_t*C_p*ρ*D_fl/(α_L*U_fl) # s, relaxation time

                T = (T_old/(α_t*Δt) + u*T_prev/Δx_fl + T_a/τ)/(1/(α_t*Δt) + u/Δx_fl + 1/τ) # °F

                T_fl_mat[n, i] = T # °F, flowline temperature
            end

            t_fl += time() # s, flowline calculation timer

            T_us = T_fl_mat[n, N_fl] # °F, upstream temperature

            !isnothing(dump_file) && @printf dump_file " %8.2f" T_us

            # end prematurely if convergence is sub-performant
            if sub_perf
                !isnothing(dump_file) && @printf dump_file " SUB PERFORMANT TIME STEP"
                cnvg = true
                break
            end

            #* choke
            t_ch -= time() # s, choke calculation timer

            p̄ = 0.5*(p_wh + p_nw) # psia, average pressure
            B = fvf(p̄, T_us)      # RB/STB, FVF
            ρ = ρ_sc/B            # lb/ft3, density
            γ = ρ/ρ_wsc           # -, specific gravity

            # q_ch = α_V*C_V*sqrt((p_wh - p_nw)/γ) # BPD, choke rate
            q_ch = α_V*C_V/B*sqrt((p_wh - p_nw)/γ) # BPD, choke rate #!!!

            t_ch += time() # s, choke calculation timer

            !isnothing(dump_file) && @printf dump_file " %8.2f" q_ch

            #* newton-raphson
            f     = q - q_ch                     # BPD, objective function
            ∂q_ch = 0.5*q_ch*∂p_wh/(p_wh - p_nw) # -, choke rate gradient
            ∂f    = 1 - ∂q_ch                    # -, objective function gradient
            λ     = 100/(j + 99)                 # -, step size
            q    -= λ*f/∂f                       # BPD, rate

            !isnothing(dump_file) && @printf dump_file " %8.1e %8.1e %8.1e %8.2f" f ∂f λ q

            # convergence check
            cnvg =
                abs(f) < ε_q && abs(p_wf - p_wf_prior) < ε_p &&
                abs(p_wh - p_wh_prior) < ε_p && abs(T_wh - T_wh_prior) < ε_T &&
                abs(T_us - T_us_prior) < ε_T
            
            !isnothing(dump_file) && @printf dump_file " %8s" cnvg
        end

        # time step data
        data.time[n]                      = t    # hr
        data.choke[n]                     = ω    # -
        data.rate[n]                      = q    # BPD
        data.bottomhole.pressure[n]       = p_wf # psia
        data.bottomhole.temperature[n]    = T_r  # °F
        data.wellhead.pressure[n]         = p_wh # psia
        data.wellhead.temperature[n]      = T_wh # °F
        data.upstream.pressure[n]         = p_wh # psia
        data.upstream.temperature[n]      = T_us # °F
        data.ambient.temperature[n]       = T_a  # °F
        data.auxilliary.iterations[n]     = j    # -
        data.auxilliary.convergence[n]    = cnvg
        data.auxilliary.subperformance[n] = sub_perf
    end

    t_sim   += time() # s, simulation timer
    t_total += time() # s, total timer

    # save timers
    data.auxilliary.timer.total          = t_total # s
    data.auxilliary.timer.setup          = t_setup # s
    data.auxilliary.timer.initialization = t_init  # s
    data.auxilliary.timer.simulation     = t_sim   # s
    data.auxilliary.timer.reservoir      = t_rsrv  # s
    data.auxilliary.timer.wellbore       = t_wb    # s
    data.auxilliary.timer.flowline       = t_fl    # s
    data.auxilliary.timer.choke          = t_ch    # s

    return data
end

