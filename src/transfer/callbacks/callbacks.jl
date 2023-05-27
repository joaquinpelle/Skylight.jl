include("syntheticpolarcap.jl")
include("black_hole_accretion_disk.jl")
include("regular_compact_object_accretion_disk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("opacities.jl")

function get_callback_and_params(configurations; kwargs...)        
    cb = get_callback(configurations.spacetime, configurations.radiative_model, coordinates_topology(configurations.spacetime))
    cb_params = get_cb_params(configurations.spacetime, configurations.radiative_model, configurations; kwargs...)
    return cb, cb_params
end
