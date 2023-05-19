#Required
set_emitter_four_velocity!(vector, position, metric, spacetime, ::AbstractRadiativeModel, coords_top) = error("set_emitter_four_velocity! not defined for this model.")
get_emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, ::AbstractRadiativeModel, coords_top) = error("get_emitted_bolometric_intensity for this model.")
get_emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, ::AbstractRadiativeModel, coords_top) = error("get_emitted_specific_intensity not defined for this model.")
is_final_position_at_source(position, spacetime, ::AbstractRadiativeModel) = error("is_final_position_at_source not defined for this model.")

#Optional
set_surface_differential!(differential, position, metric, spacetime, ::AbstractSurfaceEmissionModel, coords_top) = error("Surface differential not defined for this model.")


include("utils/utils.jl")
include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("blackholeaccretiondisk.jl")
include("bosonstaraccretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")