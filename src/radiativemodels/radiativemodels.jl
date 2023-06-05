#Required
emitter_four_velocity!(vector, position, metric, spacetime, model, coords_top) = error("emitter_four_velocity! not defined for this model.")
emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model, coords_top) = error("emitted_bolometric_intensity for this model.")
emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model, coords_top) = error("emitted_specific_intensity not defined for this model.")
is_final_position_at_source(position, spacetime, model) = error("is_final_position_at_source not defined for this model.")

#Optional for surface emission models
surface_differential!(differential, position, metric, spacetime, model, coords_top) = error("Surface differential not defined for this model.")

#For accretion disks
temperature(position, spacetime, model) = error("Temperature not defined for this model.")

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("accretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")

function unit_surface_normal!(vector, position, metric, metric_inverse, model, coords_top)
    surface_differential!(vector, position, model, coords_top)
    vector .= raise_index(vector,metric_inverse)
    normalize_spacelike!(vector, metric)
    return nothing
end