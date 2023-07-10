function perturb_data!(
    data;
    gauge_parameters=GAUGE_PARAMETERS, gauges=rand_gauges(gauge_parameters),
    corrupt_freq_fact=CORRUPT_FREQ_FACT, miss=true,
    _...
)
    # constants
    α_t = TIME_CONV_FACT # s/hr, time conversion factor

    # gauges parameters

    S1_p = gauges.pressure.short.relaxation     # s
    S2_p = gauges.pressure.short.interference   # s
    L1_p = gauges.pressure.long.relaxation      # s
    L2_p = gauges.pressure.long.initial         # psi
    M1_p = gauges.pressure.thermal.relaxation   # s
    M2_p = gauges.pressure.thermal.interference # psi-s/°F
    M3_p = gauges.pressure.thermal.deviation    # psi/°F
    N_p  = gauges.pressure.noise.deviation      # psi

    S1_T = gauges.temperature.short.relaxation   # s
    S2_T = gauges.temperature.short.interference # s
    L1_T = gauges.temperature.long.relaxation    # s
    L2_T = gauges.temperature.long.initial       # °F
    N_T  = gauges.temperature.noise.deviation    # °F

    S1_q = gauges.rate.short.relaxation     # s
    S2_q = gauges.rate.short.interference   # s
    L1_q = gauges.rate.long.relaxation      # s
    L2_q = gauges.rate.long.initial         # BPD
    M1_q = gauges.rate.thermal.relaxation   # s
    M2_q = gauges.rate.thermal.interference # BPD-s/°F
    M3_q = gauges.rate.thermal.deviation    # BPD/°F
    N_q  = gauges.rate.noise.deviation      # BPD

    # initialization

    t_prev = data.time[begin]                 # hr, previous time
    p_prev = data.upstream.pressure[begin]    # psia, previous pressure
    T_prev = data.upstream.temperature[begin] # °F, previous temperature
    q_prev = data.rate[begin]                 # BPD, previous rate

    Δp_S = 0.0  # psi, pressure short-term drift
    Δp_L = L2_p # psi, pressure long-term drift
    Δp_T = 0.0  # psi, pressure thermal drift

    ΔT_S = 0.0  # °F, temperature short-term drift
    ΔT_L = L2_T # °F, temperature long-term drift

    Δq_S = 0.0  # BPD, rate short-term drift
    Δq_L = L2_q # BPD, rate long-term drift
    Δq_T = 0.0  # BPD, rate thermal drift

    # perturbation

    data.upstream.pressure[begin]    += Δp_L # psia, pressure
    data.upstream.temperature[begin] += ΔT_L # °F, temperature
    data.rate[begin]         += Δq_L # BPD, rate

    for n in drop(eachindex(data.upstream.pressure), 1)
        t  = data.time[n]                 # hr, time
        p  = data.upstream.pressure[n]    # psia, pressure
        T  = data.upstream.temperature[n] # °F, temperature
        q  = data.rate[n]                 # BPD, rate
        Δt = α_t*(t - t_prev)             # s, time step
        
        Δp_S = (S1_p*Δp_S - S2_p*(p - p_prev))/(S1_p + Δt) # psi, pressure short-term drift
        Δp_L = Δp_L + Δt/L1_p*p                            # psi, pressure long-term drift
        Δp_T = (M1_p*Δp_T + M2_p*(T - T_prev) + M3_p*Δt*T)/(M1_p + Δt) # psi, pressure thermal drift
        Δp_N = N_p*randn()                                 # psi, pressure noise

        data.upstream.pressure[n] += Δp_S + Δp_L + Δp_T + Δp_N # psia, pressure

        ΔT_S = (S1_T*ΔT_S - S2_T*(T - T_prev))/(S1_T + Δt) # °F, temperature short-term drift
        ΔT_L = ΔT_L + Δt/L1_T*T                            # °F, temperature long-term drift
        ΔT_N = N_T*randn()                                 # °F, temperature noise
        
        data.upstream.temperature[n] += ΔT_S + ΔT_L + ΔT_N # °F, temperature
        
        Δq_S = (S1_q*Δq_S - S2_q*(q - q_prev))/(S1_q + Δt) # BPD, rate short-term drift
        Δq_L = Δq_L + Δt/L1_q*q                            # BPD, rate long-term drift
        Δq_T = (M1_q*Δq_T + M2_q*(T - T_prev) + M3_q*Δt*T)/(M1_q + Δt) # BPD, rate thermal drift
        Δq_N = N_q*randn()                                 # BPD, rate noise

        data.rate[n] += Δq_S + Δq_L + Δq_T + Δq_N # BPD, rate
        
        t_prev = t # hr, time
        p_prev = p # psia, pressure
        T_prev = T # °F, temperature
        q_prev = q # BPD, rate
    end

    miss || return data

    # missing

    permutations = (
        # three readings
        (true, true, true),
        (true, true, true),
        (true, true, true),
        # two readings
        (true, true, false),
        (true, false, true),
        (false, true, true),
        # one reading
        (true, false, false),
        (false, true, false),
        (false, false, true)
    )

    n = length(data.upstream.pressure)
    idx_last_change = firstindex(data.rate)
    pres_valid, temp_valid, rate_valid = rand(permutations)

    for i in eachindex(data.rate)
        change_prob = ((i - idx_last_change)/n)^corrupt_freq_fact

        if rand() < change_prob
            pres_valid, temp_valid, rate_valid = rand(permutations)
            idx_last_change = i
        end

        data.upstream.pressure[i]    = ifelse(pres_valid, data.upstream.pressure[i], NaN)
        data.upstream.temperature[i] = ifelse(temp_valid, data.upstream.temperature[i], NaN)
        data.rate[i]         = ifelse(rate_valid, data.rate[i], NaN)
    end

    return data
end
