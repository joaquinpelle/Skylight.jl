export get_callback_and_params

include("callbacks.jl")
include("cb_params.jl")

function get_callback_and_params(configurations)
        
    cb = get_callback(configurations.emission_model, coordinate_system_class(configurations.emission_model))
    cb_params = get_cb_params(configurations.emission_model, configurations)

    return cb, cb_params

end