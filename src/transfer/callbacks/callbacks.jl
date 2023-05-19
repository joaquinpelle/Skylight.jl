include("syntheticpolarcap.jl")
include("blackholeaccretiondisk.jl")
include("bosonstaraccretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("opacities.jl")


function get_callback_and_params(configurations; kwargs...)
        
    cb = get_callback(configurations.radiative_model, coordinates_topology(configurations.spacetime))
    cb_params = get_cb_params(configurations.radiative_model, configurations; kwargs...)

    return cb, cb_params

end

get_rmax(configurations::AbstractETOConfigurations) = configurations.observer_distance

function get_rmax(configurations::AbstractOTEConfigurations) 
    
    d = configurations.image_plane.distance
    hs = configurations.image_plane.horizontal_side_image_plane
    vs = configurations.image_plane.vertical_side_image_plane
    return 1.1*sqrt(d^2 + vs^2 + hs^2)

end