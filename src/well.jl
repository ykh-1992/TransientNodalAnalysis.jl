struct Oil
    gravity::Float64         # °API
    viscosity::Float64       # cP
    compressibility::Float64 # 1/psi
    expansion::Float64       # 1/°F
    capacity::Float64        # Btu/lb-°F
    # factor::Float64          # RB/STB #!!!
end

struct Reservoir
    pressure::Float64        # psia
    permeability::Float64    # mD
    porosity::Float64        # -
    thickness::Float64       # ft
    radius::Float64          # ft
    compressibility::Float64 # 1/psi
    temperature::Float64     # °F
    skin::Float64            # -
end

abstract type Boundary end

function transient_function(::Boundary, H, r_w, t) end

struct Wellbore
    length::Float64    # ft
    diameter::Float64  # in
    transfer::Float64  # Btu/hr-ft2-°F
    roughness::Float64 # ft
end

struct Flowline
    length::Float64   # ft
    diameter::Float64 # in
    transfer::Float64 # Btu/hr-ft2-°F
end

struct Choke
    coefficient::Float64 # gpm/sqrt(psi), maximum flow coefficient
    curvature::Float64   # -, curvature parameter
    inflection::Float64  # -, inflection point parameter
end

function flow_coefficient(choke, opening)
    m = choke.coefficient # gpm/sqrt(psi), maximum flow coefficient
    a = choke.curvature   # -, curvature parameter
    b = choke.inflection  # -, inflection point parameter
    ω = opening           # -, choke opening

    s = 1/(1/(1 + exp(a*(b - 1))) - 1/(1 + exp(a*b))) # -, scaling parameter
    c = 1/(1 + exp(a*b))                              # -, intercept parameter
    y = m*s*(1/(1 + exp(a*(b - ω))) - c)              # gpm/sqrt(psi), flow coefficient

    return y
end

struct Ambient
    average::Float64 # °F
    annual::Float64  # °F
    daily::Float64   # °F
    shift::Float64   # -
end

struct Network
    pressure::Float64 # psia
end

struct Well{B<:Boundary}
    oil::Oil
    reservoir::Reservoir
    boundary::B
    wellbore::Wellbore
    flowline::Flowline
    choke::Choke
    ambient::Ambient
    network::Network
end
