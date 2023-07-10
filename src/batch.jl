using Printf
using OffsetArrays: Origin
using DataFrames
using CSV
using .TransientNodalAnalysis

const WELLS_HEADER = [
    "well",
    "oil.gravity",
    "oil.viscosity",
    "oil.compressibility",
    "oil.expansion",
    "oil.capacity",
    "reservoir.pressure",
    "reservoir.permeability",
    "reservoir.porosity",
    "reservoir.thickness",
    "reservoir.radius",
    "reservoir.compressibility",
    "reservoir.temperature",
    "reservoir.skin",
    "boundary.type",
    "boundary.distance",
    "wellbore.length",
    "wellbore.diameter",
    "wellbore.transfer",
    "wellbore.roughness",
    "flowline.length",
    "flowline.diameter",
    "flowline.transfer",
    "choke.coefficient",
    "choke.curvature",
    "choke.inflection",
    "ambient.average",
    "ambient.annual",
    "ambient.daily",
    "ambient.shift",
    "network.pressure"
]

const AUX_HEADER = [
    "well",
    "iterations.total",
    "iterations.max",
    "steps.nonconvergence",
    "steps.subperformance",
    "timer.total",
    "timer.setup",
    "timer.initialization",
    "timer.simulation",
    "timer.simulation.reservoir",
    "timer.simulation.wellbore",
    "timer.simulation.flowline",
    "timer.simulation.choke"
]

const DATA_HEADER = [
    "well",
    "step",
    "time",
    "choke",
    "rate",
    "flowing",
    "bottomhole.pressure",
    "bottomhole.temperature",
    "wellhead.pressure",
    "wellhead.temperature",
    "wellbore.friction",
    "upstream.pressure",
    "upstream.temperature",
    "ambient.temperature",
    "convergence",
    "iterations",
    "subperformance"
]

function run_batch_simulation(
    num_wells;
    num_time_step=NUM_TIME_STEP,
    parameters=WELL_PARAMETERS,
    max_wb_seg_len=MAX_WB_SEG_LEN,
    max_fl_seg_len=MAX_FL_SEG_LEN,
    path=".",
    kwargs...
)
    t_all = -time()

    well_vec   = Vector{Well}(undef, num_wells)
    max_wb_len = 0.0
    max_fl_len = 0.0

    for i in eachindex(well_vec)
        well        = rand_well(parameters)
        max_wb_len  = max(max_wb_len, well.wellbore.length)
        max_fl_len  = max(max_fl_len, well.flowline.length)
        well_vec[i] = well
    end

    num_wb_seg = ceil(Int, max_wb_len/max_wb_seg_len)
    num_fl_seg = ceil(Int, max_fl_len/max_fl_seg_len)

    num_threads = Threads.nthreads()

    choke_vec = [zeros(num_time_step) for _ in 1:num_threads]
    data_vec  = [SimulationData(num_time_step, num_wb_seg, num_fl_seg) for _ in 1:num_threads]

    wells_dataframe = DataFrame(
        "well"                      => 1:num_wells,
        "oil_gravity"               => [well.oil.gravity for well in well_vec],
        "oil_viscosity"             => [well.oil.viscosity for well in well_vec],
        "oil_compressibility"       => [well.oil.compressibility for well in well_vec],
        "oil_expansion"             => [well.oil.expansion for well in well_vec],
        "oil_capacity"              => [well.oil.capacity for well in well_vec],
        "reservoir_pressure"        => [well.reservoir.pressure for well in well_vec],
        "reservoir_permeability"    => [well.reservoir.permeability for well in well_vec],
        "reservoir_porosity"        => [well.reservoir.porosity for well in well_vec],
        "reservoir_thickness"       => [well.reservoir.thickness for well in well_vec],
        "reservoir_radius"          => [well.reservoir.radius for well in well_vec],
        "reservoir_compressibility" => [well.reservoir.compressibility for well in well_vec],
        "reservoir_temperature"     => [well.reservoir.temperature for well in well_vec],
        "reservoir_skin"            => [well.reservoir.skin for well in well_vec],
        "boundary_type"             => [boundary_name(well.boundary) for well in well_vec],
        "boundary_distance"         => [boundary_distance(well.boundary) for well in well_vec],
        "wellbore_length"           => [well.wellbore.length for well in well_vec],
        "wellbore_diameter"         => [well.wellbore.diameter for well in well_vec],
        "wellbore_transfer"         => [well.wellbore.transfer for well in well_vec],
        "wellbore_roughness"        => [well.wellbore.roughness for well in well_vec],
        "flowline_length"           => [well.flowline.length for well in well_vec],
        "flowline_diameter"         => [well.flowline.diameter for well in well_vec],
        "flowline_transfer"         => [well.flowline.transfer for well in well_vec],
        "choke_coefficient"         => [well.choke.coefficient for well in well_vec],
        "choke_curvature"           => [well.choke.curvature for well in well_vec],
        "choke_inflection"          => [well.choke.inflection for well in well_vec],
        "ambient_average"           => [well.ambient.average for well in well_vec],
        "ambient_annual"            => [well.ambient.annual for well in well_vec],
        "ambient_daily"             => [well.ambient.daily for well in well_vec],
        "ambient_shift"             => [well.ambient.shift for well in well_vec],
        "network_pressure"          => [well.network.pressure for well in well_vec]
    )

    num_rows = num_wells*(num_time_step + 1)

    data_dataframe = DataFrame(
        "well"                   => Vector{Int}(undef, num_rows),
        "step"                   => Vector{Int}(undef, num_rows),
        "time"                   => Vector{Float64}(undef, num_rows),
        "choke"                  => Vector{Float64}(undef, num_rows),
        "rate"                   => Vector{Float64}(undef, num_rows),
        "flowing"                => BitVector(undef, num_rows),
        "bottomhole_pressure"    => Vector{Float64}(undef, num_rows),
        "bottomhole_temperature" => Vector{Float64}(undef, num_rows),
        "wellhead_pressure"      => Vector{Float64}(undef, num_rows),
        "wellhead_temperature"   => Vector{Float64}(undef, num_rows),
        "wellbore_friction"      => Vector{Float64}(undef, num_rows),
        "upstream_pressure"      => Vector{Float64}(undef, num_rows),
        "upstream_temperature"   => Vector{Float64}(undef, num_rows),
        "ambient_temperature"    => Vector{Float64}(undef, num_rows),
        "convergence"            => Vector{Float64}(undef, num_rows),
        "iterations"             => Vector{Float64}(undef, num_rows),
        "subperformance"         => Vector{Float64}(undef, num_rows)
    )

    aux_dataframe = DataFrame(
        "well"                       => 1:num_wells,
        "iterations_total"           => Vector{Int}(undef, num_wells),
        "iterations_max"             => Vector{Int}(undef, num_wells),
        "steps_nonconvergence"       => Vector{Int}(undef, num_wells),
        "steps_subperformance"       => Vector{Int}(undef, num_wells),
        "timer_total"                => Vector{Float64}(undef, num_wells),
        "timer_setup"                => Vector{Float64}(undef, num_wells),
        "timer_initialization"       => Vector{Float64}(undef, num_wells),
        "timer_simulation"           => Vector{Float64}(undef, num_wells),
        "timer_simulation_reservoir" => Vector{Float64}(undef, num_wells),
        "timer_simulation_wellbore"  => Vector{Float64}(undef, num_wells),
        "timer_simulation_flowline"  => Vector{Float64}(undef, num_wells),
        "timer_simulation_choke"     => Vector{Float64}(undef, num_wells)
    )

    t_wall    = -time()
    t_threads = [0.0 for _ in 1:num_threads]

    Threads.@threads for i in eachindex(well_vec)
        k = Threads.threadid()

        t_threads[k] -= time()

        choke = choke_vec[k]
        data  = data_vec[k]
        well  = well_vec[i]

        rand_choke_opening!(choke; kwargs...)
        simulate_well!(
            data, well, choke;
            max_wb_seg_len, max_fl_seg_len, kwargs...
        )

        for n in eachindex(data.time)
            j = (i - 1)*(num_time_step + 1) + n + 1
            data_dataframe.flowing[j] = data.rate[n] != 0
        end

        perturb_data!(data)

        iterations      = 0
        max_iterations  = 0
        non_convergence = 0
        subperformance  = 0

        for n in eachindex(data.time)
            iterations      += data.auxilliary.iterations[n]
            max_iterations   = max(data.auxilliary.iterations[n], max_iterations)
            non_convergence += !data.auxilliary.convergence[n]
            subperformance  += data.auxilliary.subperformance[n]

            j = (i - 1)*(num_time_step + 1) + n + 1

            friction = sum(data.wellbore.friction[n, m] for m in axes(data.wellbore.friction, 2))

            data_dataframe.well[j]                   = i
            data_dataframe.step[j]                   = n
            data_dataframe.time[j]                   = data.time[n]
            data_dataframe.choke[j]                  = data.choke[n]
            data_dataframe.rate[j]                   = data.rate[n]
            data_dataframe.bottomhole_pressure[j]    = data.bottomhole.pressure[n]
            data_dataframe.bottomhole_temperature[j] = data.bottomhole.temperature[n]
            data_dataframe.wellhead_pressure[j]      = data.wellhead.pressure[n]
            data_dataframe.wellhead_temperature[j]   = data.wellhead.temperature[n]
            data_dataframe.wellbore_friction[j]      = friction
            data_dataframe.upstream_pressure[j]      = data.upstream.pressure[n]
            data_dataframe.upstream_temperature[j]   = data.upstream.temperature[n]
            data_dataframe.ambient_temperature[j]    = data.ambient.temperature[n]
            data_dataframe.convergence[j]            = data.auxilliary.convergence[n]
            data_dataframe.iterations[j]             = data.auxilliary.iterations[n]
            data_dataframe.subperformance[j]         = data.auxilliary.subperformance[n]
        end

        aux_dataframe.well[i]                       = i
        aux_dataframe.iterations_total[i]           = iterations
        aux_dataframe.iterations_max[i]             = max_iterations
        aux_dataframe.steps_nonconvergence[i]       = non_convergence
        aux_dataframe.steps_subperformance[i]       = subperformance
        aux_dataframe.timer_total[i]                = data.auxilliary.timer.total
        aux_dataframe.timer_setup[i]                = data.auxilliary.timer.setup
        aux_dataframe.timer_initialization[i]       = data.auxilliary.timer.initialization
        aux_dataframe.timer_simulation[i]           = data.auxilliary.timer.simulation
        aux_dataframe.timer_simulation_reservoir[i] = data.auxilliary.timer.reservoir
        aux_dataframe.timer_simulation_wellbore[i]  = data.auxilliary.timer.wellbore
        aux_dataframe.timer_simulation_flowline[i]  = data.auxilliary.timer.flowline
        aux_dataframe.timer_simulation_choke[i]     = data.auxilliary.timer.choke

        @printf "Thread %d completed well %04d (%.3f sec simulation time)\n" k i data.auxilliary.timer.total
        
        t_threads[k] += time()
    end

    println("\nCompleted simulation")
    
    t_wall += time()
    @printf "CPU time: %.3f sec\n" sum(t_threads)
    @printf "Wall time: %.3f sec\n\n" t_wall

    t_csv = -time()
    
    wells_path = joinpath(path, "wells.csv")
    data_path  = joinpath(path, "data.csv")
    aux_path   = joinpath(path, "auxilliary.csv")

    CSV.write(wells_path, wells_dataframe)
    CSV.write(aux_path, aux_dataframe)
    CSV.write(data_path, data_dataframe)

    t_csv += time()
    println("Completed CSV writing")
    @printf "CSV time: %.3f sec\n\n" t_csv

    t_all += time()
    @printf "Total time: %.3f sec\n" t_all
end

joincomma(X...) = join(X, ',')

boundary_name(boundary::IARF)                     = "IARF"
boundary_name(boundary::NoFlow)                   = "NoFlow"
boundary_name(boundary::ParallelNoFlow)           = "ParallelNoFlow"
boundary_name(boundary::BoxNoFlow)                = "BoxNoFlow"
boundary_name(boundary::ConstantPressure)         = "ConstantPressure"
boundary_name(boundary::ParallelConstantPressure) = "ParallelConstantPressure"
boundary_name(boundary::BoxConstantPressure)      = "BoxConstantPressure"

boundary_distance(boundary)                              = boundary.distance
boundary_distance(boundary::IARF) = 0.0

run_batch_simulation(100)
