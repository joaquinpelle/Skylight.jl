#Required
set_emitter_four_velocity!(vector, position, metric, spacetime, ::AbstractRadiativeModel, coords_top) = error("set_emitter_four_velocity! not defined for this model.")
get_emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, ::AbstractRadiativeModel, coords_top) = error("get_emitted_bolometric_intensity for this model.")
get_emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, ::AbstractRadiativeModel, coords_top) = error("get_emitted_specific_intensity not defined for this model.")
is_final_position_at_source(position, spacetime, ::AbstractRadiativeModel) = error("is_final_position_at_source not defined for this model.")

#Optional
set_surface_differential!(differential, position, metric, spacetime, ::AbstractSurfaceEmissionModel, coords_top) = error("Surface differential not defined for this model.")

#For accretion disks
temperature(position, spacetime, ::AbstractAccretionDisk) = error("Temperature not defined for this model.")

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("accretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")

function set_unit_surface_normal!(vector, position, metric, metric_inverse, model, coords_top)
    set_surface_differential!(vector, position, model, coords_top)
    vector .= raise_index(vector,metric_inverse)
    normalize_spacelike!(vector, metric)
end