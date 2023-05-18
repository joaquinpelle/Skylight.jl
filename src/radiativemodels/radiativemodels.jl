#Required
set_emitter_four_velocity!(vector, position, metric, spacetime, model::AbstractRadiativeModel, coord_system) = error("set_emitter_four_velocity! not defined for this model.")
get_emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model::AbstractRadiativeModel, coord_system) = error("get_emitted_bolometric_intensity for this model.")
get_emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model::AbstractRadiativeModel, coord_system) = error("get_emitted_specific_intensity not defined for this model.")
is_final_position_at_source(position, spacetime, model::AbstractRadiativeModel) = error("is_final_position_at_source not defined for this model.")

#Optional
set_surface_differential!(differential, position, metric, spacetime, model::SurfaceEmissionModel, coord_system) = error("Surface differential not defined for this model.")


include("utils/utils.jl")
include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("blackholeaccretiondisk.jl")
include("bosonstaraccretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")

#Required
export set_emitter_four_velocity!
export get_emitted_bolometric_intensity
export get_emitted_specific_intensity
export is_final_position_at_source
export set_surface_differential!
