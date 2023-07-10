function rand_well(params=WELL_PARAMETERS; well_parameters=params, _...)
    return Well(
        Oil(
            rand(well_parameters.oil.gravity),         # °API
            rand(well_parameters.oil.viscosity),       # cP
            rand(well_parameters.oil.compressibility), # 1/psi
            rand(well_parameters.oil.expansion),       # 1/°F
            rand(well_parameters.oil.capacity)         # Btu/lb-°F
        ),
        Reservoir(
            rand(well_parameters.reservoir.pressure),        # psia
            rand(well_parameters.reservoir.permeability),    # mD
            rand(well_parameters.reservoir.porosity),        # -
            rand(well_parameters.reservoir.thickness),       # ft
            rand(well_parameters.reservoir.radius),          # ft
            rand(well_parameters.reservoir.compressibility), # 1/psi
            rand(well_parameters.reservoir.temperature),     # °F
            rand(well_parameters.reservoir.skin)             # -
        ),
        rand(well_parameters.boundary.types)(
            rand(well_parameters.boundary.distance) # ft
        ),
        Wellbore(
            rand(well_parameters.wellbore.length),   # ft
            rand(well_parameters.wellbore.diameter), # in
            rand(well_parameters.wellbore.transfer), # Btu/hr-ft2-°F
            rand(well_parameters.wellbore.roughness) # ft
        ),
        Flowline(
            rand(well_parameters.flowline.length),   # ft
            rand(well_parameters.flowline.diameter), # in
            rand(well_parameters.flowline.transfer)  # Btu/hr-ft2-°F
        ),
        Choke(
            rand(well_parameters.choke.coefficient), # gpm/sqrt(psi)
            rand(well_parameters.choke.curvature),   # -
            rand(well_parameters.choke.inflection)   # -
        ),
        Ambient(
            rand(well_parameters.ambient.average), # °F
            rand(well_parameters.ambient.annual),  # °F
            rand(well_parameters.ambient.daily),   # °F
            rand(well_parameters.ambient.shift)    # -
        ),
        Network(
            rand(well_parameters.network.pressure) # psia
        )
    )
end

function rand_choke_opening(
    n=NUM_TIME_STEP;
    num_time_step=n,
    shutin_prob=SHUTIN_PROB,
    choke_freq_fact=CHOKE_FREQ_FACT,
    _...
)
    choke_opening = zeros(num_time_step)
    return rand_choke_opening!(choke_opening; shutin_prob, choke_freq_fact)
end

function rand_choke_opening!(
    choke_opening;
    shutin_prob=SHUTIN_PROB,
    choke_freq_fact=CHOKE_FREQ_FACT,
    _...
)
    rand_opening() = rand() < shutin_prob ? 0.0 : rand()

    n = length(choke_opening)     # -, number of steps
    ω = rand_opening()            # -, choke opening
    j = firstindex(choke_opening) # -, index of last choke opening

    for i in eachindex(choke_opening)
        change_prob = ((i - j)/n)^choke_freq_fact # -, probability to change choke opening

        if rand() < change_prob
            ω = rand_opening() # -, choke opening
            j = i              # -, index of last choke opening
        end

        choke_opening[i] = ω # -, choke opening
    end

    return choke_opening
end

# TODO: generate gauges as structs to be saved and used later
function rand_gauges(params=GAUGE_PARAMETERS; gauge_parameters=params, _...)
    return (
        pressure = (
            short = (
                relaxation   = rand(gauge_parameters.pressure.short.relaxation),  # s
                interference = rand(gauge_parameters.pressure.short.interference) # s
            ),
            long = (
                relaxation = rand(gauge_parameters.pressure.long.relaxation), # s
                initial    = rand(gauge_parameters.pressure.long.initial)     # s
            ),
            thermal = (
                relaxation   = rand(gauge_parameters.pressure.thermal.relaxation),   # s
                interference = rand(gauge_parameters.pressure.thermal.interference), # psi-s/°F
                deviation    = rand(gauge_parameters.pressure.thermal.deviation)     # psi/°F
            ),
            noise = (
                deviation = rand(gauge_parameters.pressure.noise.deviation), # psi
            )
        ),
        temperature = (
            short = (
                relaxation   = rand(gauge_parameters.temperature.short.relaxation),  # s
                interference = rand(gauge_parameters.temperature.short.interference) # s
            ),
            long = (
                relaxation = rand(gauge_parameters.temperature.long.relaxation), # s
                initial    = rand(gauge_parameters.temperature.long.initial)     # s
            ),
            noise = (
                deviation = rand(gauge_parameters.temperature.noise.deviation), # °F
            )
        ),
        rate = (
            short = (
                relaxation   = rand(gauge_parameters.rate.short.relaxation),  # s
                interference = rand(gauge_parameters.rate.short.interference) # s
            ),
            long = (
                relaxation = rand(gauge_parameters.rate.long.relaxation), # s
                initial    = rand(gauge_parameters.rate.long.initial)     # s
            ),
            thermal = (
                relaxation   = rand(gauge_parameters.rate.thermal.relaxation),   # s
                interference = rand(gauge_parameters.rate.thermal.interference), # BPD-s/°F
                deviation    = rand(gauge_parameters.rate.thermal.deviation)     # BPD/°F
            ),
            noise = (
                deviation = rand(gauge_parameters.rate.noise.deviation), # BPD
            )
        )
    )
end
