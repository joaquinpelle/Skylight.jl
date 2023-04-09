export get_callback_and_params

abstract type CallbackParameters end

include("neutronstarhotspots.jl")
include("blackholeaccretiondisk.jl")
include("bosonstaraccretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("opacities.jl")


function get_callback_and_params(configurations; kwargs...)
        
    cb = get_callback(configurations.radiative_model, coordinate_system_class(configurations.spacetime))
    cb_params = get_cb_params(configurations.radiative_model, configurations; kwargs...)

    return cb, cb_params

end

get_rmax(configurations::ETOConfigurations) = configurations.observer_distance

function get_rmax(configurations::OTEConfigurations) 
    
    d = configurations.image_plane.distance
    hs = configurations.image_plane.horizontal_side_image_plane
    vs = configurations.image_plane.vertical_side_image_plane
    return 1.1*sqrt(d^2 + vs^2 + hs^2)

end