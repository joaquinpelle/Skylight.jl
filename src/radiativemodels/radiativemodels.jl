#Required
emitter_four_velocity!(vector, position, metric, spacetime, model, coords_top, cache) = error("emitter_four_velocity! not defined for this model.")
emitter_four_velocity!(vector, position, metric, spacetime, model, coords_top) = error("emitter_four_velocity! not defined for this model.")
emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model, coords_top, cache) = error("emitted_bolometric_intensity for this model.")
emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model, coords_top) = error("emitted_bolometric_intensity for this model.")
emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model, coords_top, cache) = error("emitted_specific_intensity not defined for this model.")
emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model, coords_top) = error("emitted_specific_intensity not defined for this model.")
line_emission_profile(position, momentum, emitter_four_velocity, metric, spacetime, model, coords_top, cache) = error("line_emission_profile not defined for this model.")
is_final_position_at_source(position, spacetime, model) = error("is_final_position_at_source not defined for this model.")

#Optional for surface emission models
allocate_cache(::AbstractRadiativeModel) = nothing
surface_differential!(differential, position, model, coords_top) = error("Surface differential not defined for this model.")

#For accretion disks
temperature(position, spacetime, model) = error("Temperature not defined for this model.")

include("radiativeprocesses/thermalemission.jl")
include("radiativeprocesses/bremsstrahlung.jl")
include("radiativeprocesses/synchrotron.jl")

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("accretiondisk.jl")
include("staracrosswormhole.jl")
include("verticalscreen.jl")
include("lamppostcorona.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")

function unit_surface_normal!(vector, position, metric, metric_inverse, model, coords_top)
    surface_differential!(vector, position, model, coords_top)
    vector .= raise_index(vector,metric_inverse)
    normalize_spacelike!(vector, metric)
    return nothing
end

emitter_four_velocity!(v, position, metric, spacetime, model, coords_top, ::Nothing) = emitter_four_velocity!(v, position, metric, spacetime, model, coords_top)
emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model, coords_top, ::Nothing) = emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model, coords_top)
emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model, coords_top, ::Nothing) = emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model, coords_top)