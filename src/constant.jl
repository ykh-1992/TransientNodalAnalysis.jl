const STD_COND_PRES     = 14.7  # psia, standard condition pressure
const STD_COND_TEMP     = 60    # °F, standard condition temperature
const STD_COND_WTR_DENS = 62.4  # lb/ft3, standard condition water density
const GRAV_ACC          = 32.17 # ft/s2, gravitational acceleration

const LEN_CONV_FACT   = 12         # in/ft, length conversion factor
const TIME_CONV_FACT  = 3600       # s/hr, time conversion factor
const GRAV_CONV_FACT  = 32.17      # lbm-ft/lbf-s2, gravity conversion factor
const RATE_CONV_FACT  = 15388.5    # BPD/ft3/s, flow rate conversion factor
const VISC_CONV_FACT  = 6894757    # cP/psi-s, viscosity conversion factor
const PERM_CONV_FACT  = 9.41363e13 # mD/ft2, permeability conversion factor
const VALVE_CONV_FACT = 34.2857    # BPD/gpm, valve flow rate conversion factor

const SHUTIN_PROB     = 0.5 # -, probability of well being shut-in
const CHOKE_FREQ_FACT = 1.4 # -, frequency factor for changing choke opening

const NUM_TIME_STEP  = 43800 # -, number of time steps
const TIME_STEP      = 0.2   # hr, time step
const MAX_ITER       = 1000  # -, maximum number of convergence iterations
const RATE_ABS_TOL   = 1e-1  # BPD, absolute tolerance on rate
const PRES_ABS_TOL   = 1e-1  # psi, absolute tolerance on pressure
const TEMP_ABS_TOL   = 1e-1  # °F, absolute tolerance on temperature
const MAX_WB_SEG_LEN = 50    # ft, maximum wellbore segment length
const MAX_FL_SEG_LEN = 3     # ft, maximum flowline segment length
const NEAR_ZERO      = 1e-16 # -, near zero value

const RATE_REL_DEV = 0.003 # -, rate relative noise deviation
const RATE_ABS_DEV = 5     # BPD, rate absolute noise deviation
const PRES_REL_DEV = 0.001 # -, pressure relative noise deviation
const PRES_ABS_DEV = 2     # psi, pressure absolute noise deviation
const TEMP_REL_DEV = 0.001 # -, temperature relative noise deviation
const TEMP_ABS_DEV = 0.1   # °F, temperature absolute noise deviation

const CORRUPT_FREQ_FACT = 3 # -, frequency factor for injecting missing or frozen data

# well parameters
const WELL_PARAMETERS = (
    oil = (
        gravity         = Uniform(40, 70),     # °API
        viscosity       = Uniform(1, 10),      # cP
        compressibility = Uniform(5e-6, 1e-4), # 1/psi
        expansion       = Uniform(1e-4, 5e-4), # 1/°F
        capacity        = Uniform(0.5, 1.0)    # Btu/lb-°F
    ),
    reservoir = (
        pressure        = Uniform(3000, 8000),     # psia
        permeability    = Uniform(10, 200),        # mD
        porosity        = Uniform(0.1, 0.3),       # -
        thickness       = Uniform(5, 100),         # ft
        radius          = Uniform(0.0833, 0.3333), # ft
        compressibility = Uniform(1e-6, 1e-5),     # 1/psi
        temperature     = Uniform(120, 250),       # °F
        skin            = Uniform(-2, 10)          # -
    ),
    boundary = (
        types = [
            IARF, 
            NoFlow, ParallelNoFlow, BoxNoFlow, 
            ConstantPressure, ParallelConstantPressure, BoxConstantPressure
        ],
        distance = Uniform(500, 10000) # ft
    ),
    wellbore = (
        length    = Uniform(2000, 8000), # ft
        diameter  = Uniform(1, 5),       # in
        transfer  = Uniform(0.2, 1.25),  # Btu/hr-ft2-°F
        roughness = Uniform(5e-5, 1e-3)  # ft
    ),
    flowline = (
        length   = Uniform(10, 200),   # ft
        diameter = Uniform(1, 3),      # in
        transfer = Uniform(0.833, 2.5) # Btu/hr-ft2-°F
    ),
    ambient = (
        average = Uniform(50, 100), # °F
        annual  = Uniform(5, 20),   # °F
        daily   = Uniform(5, 15),   # °F
        shift   = Uniform(0, 2*π)   # -
    ),
    choke = (
        coefficient = Uniform(10, 25),  # gpm/sqrt(psi)
        curvature   = Uniform(1, 10),   # -
        inflection  = Uniform(0.7, 0.8) # -
    ),
    network = (
        pressure = Uniform(100, 1000), # psia
    )
)

# gauge parameters
const GAUGE_PARAMETERS = (
    pressure = (
        short = (
            relaxation   = Uniform(0, 120), # s
            interference = Uniform(0, 6)    # s
        ),
        long = (
            relaxation = Uniform(1e8, 1e10), # s
            initial    = Uniform(0, 5)       # psi
        ),
        thermal = (
            relaxation   = Uniform(0, 410),         # s
            interference = Uniform(0, 169),         # psi-s/°F
            deviation    = Uniform(-0.0172, 0.0154) # psi/°F
        ),
        noise = (
            deviation = Uniform(0, 10), # psi
        )
    ),
    temperature = (
        short = (
            relaxation   = Uniform(0, 120), # s
            interference = Uniform(0, 6)    # s
        ),
        long = (
            relaxation = Uniform(1e9, 1e10), # s
            initial    = Uniform(0, 1)       # °F
        ),
        noise = (
            deviation = Uniform(0, 2), # °F
        )
    ),
    rate = (
        short = (
            relaxation   = Uniform(0, 120), # s
            interference = Uniform(0, 6)    # s
        ),
        long = (
            relaxation = Uniform(1e8, 1e9), # s
            initial    = Uniform(0, 10)     # BPD
        ),
        thermal = (
            relaxation   = Uniform(0, 410),         # s
            interference = Uniform(0, 169),         # BPD-s/°F
            deviation    = Uniform(-0.0172, 0.0154) # BPD/°F
        ),
        noise = (
            deviation = Uniform(0, 15), # BPD
        )
    )
)
