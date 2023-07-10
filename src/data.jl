oazeros(T::DataType, dims...) = OffsetArrays.Origin(0)(zeros(T, (dim + 1 for dim in dims)...))
oazeros(dims...)            = oazeros(Float64, dims...)

struct SurfaceData
    rate::OffsetVector{Float64, Vector{Float64}}  # BPD
    choke::OffsetVector{Float64, Vector{Float64}} # -
end

struct NodeData
    pressure::OffsetVector{Float64, Vector{Float64}}    # psia
    temperature::OffsetVector{Float64, Vector{Float64}} # °F
end

NodeData(n) = NodeData(oazeros(n), oazeros(n))

Base.copy(data::NodeData) = NodeData(copy(data.pressure), copy(data.temperature))

struct TemperatureNodeData
    temperature::OffsetVector{Float64, Vector{Float64}} # °F
end

TemperatureNodeData(n) = TemperatureNodeData(oazeros(n))

Base.copy(data::TemperatureNodeData) = TemperatureNodeData(copy(data.temperature))

struct WellboreSegmentData
    pressure::OffsetMatrix{Float64, Matrix{Float64}}    # psia
    temperature::OffsetMatrix{Float64, Matrix{Float64}} # °F
    friction::OffsetMatrix{Float64, Matrix{Float64}}    # psi
end

WellboreSegmentData(n_t, n_x) = WellboreSegmentData(
    oazeros(n_t, n_x), oazeros(n_t, n_x), oazeros(n_t, n_x)
)

Base.copy(data::WellboreSegmentData) = WellboreSegmentData(
    copy(data.pressure), copy(data.temperature), copy(data.friction)
)

struct FlowlineSegmentData
    temperature::OffsetMatrix{Float64, Matrix{Float64}} # °F
end

FlowlineSegmentData(n_t, n_x) = FlowlineSegmentData(oazeros(n_t, n_x))

Base.copy(data::FlowlineSegmentData) = FlowlineSegmentData(copy(data.temperature))

mutable struct TimerData
    total::Float64          # s
    setup::Float64          # s
    initialization::Float64 # s
    simulation::Float64     # s
    reservoir::Float64      # s
    wellbore::Float64       # s
    flowline::Float64       # s
    choke::Float64          # s
end

TimerData() = TimerData(0, 0, 0, 0, 0, 0, 0, 0)

Base.copy(timer::TimerData) = TimerData(
    timer.total, timer.setup, timer.initialization, timer.simulation,
    timer.reservoir, timer.wellbore, timer.flowline, timer.choke
)

struct AuxilliaryData
    transient::OffsetVector{Float64, Vector{Float64}} # -
    convergence::OffsetVector{Bool, Vector{Bool}}
    iterations::OffsetVector{Float64, Vector{Float64}}
    subperformance::OffsetVector{Bool, Vector{Bool}}
    timer::TimerData
end

AuxilliaryData(n) = AuxilliaryData(
    oazeros(n), oazeros(Bool, n), oazeros(Int, n), oazeros(Bool, n), TimerData()
)

Base.copy(data::AuxilliaryData) = AuxilliaryData(
    copy(data.transient), copy(data.convergence), copy(data.iterations), copy(data.subperformance),
    copy(data.timer)
)

struct SimulationData
    time::OffsetVector{Float64, Vector{Float64}}  # hr
    choke::OffsetVector{Float64, Vector{Float64}} # -
    rate::OffsetVector{Float64, Vector{Float64}}  # BPD
    bottomhole::NodeData
    wellbore::WellboreSegmentData
    wellhead::NodeData
    flowline::FlowlineSegmentData
    upstream::NodeData
    ambient::TemperatureNodeData
    auxilliary::AuxilliaryData
end

SimulationData(n_t, n_wb, n_fl) = SimulationData(
    oazeros(n_t),
    oazeros(n_t),
    oazeros(n_t),
    NodeData(n_t),
    WellboreSegmentData(n_t, n_wb),
    NodeData(n_t),
    FlowlineSegmentData(n_t, n_fl),
    NodeData(n_t),
    TemperatureNodeData(n_t),
    AuxilliaryData(n_t)
)

function SimulationData(
    well;
    num_time_step=NUM_TIME_STEP,
    max_wb_seg_len=MAX_WB_SEG_LEN, max_fl_seg_len=MAX_FL_SEG_LEN
)
    n_wb = ceil(Int, well.wellbore.length/max_wb_seg_len) # -, number of wellbore segments
    n_fl = ceil(Int, well.flowline.length/max_fl_seg_len) # -, number of flowline segments
    data = SimulationData(num_time_step, n_wb, n_fl)

    return data
end

Base.copy(data::SimulationData) = SimulationData(
    copy(data.time), copy(data.choke), copy(data.rate), copy(data.bottomhole), copy(data.wellbore),
    copy(data.wellhead), copy(data.flowline), copy(data.upstream), copy(data.ambient),
    copy(data.auxilliary)
)
