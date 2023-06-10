include("syntheticpolarcap.jl")
include("black_hole_accretion_disk.jl")
include("regular_compact_object_accretion_disk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("opacities.jl")

function callback_setup(configurations; kwargs...)        
    cb = callback(configurations.spacetime, configurations.radiative_model, coordinates_topology(configurations.spacetime))
    cbp = callback_parameters(configurations.spacetime, configurations.radiative_model, configurations; kwargs...)
    return cb, cbp
end

#Interface
callback(spacetime, model, configurations; kwargs...) = throw(ArgumentError("Callback not implemented for this spacetime and radiative model"))
callback_parameters(spacetime, model, configurations; kwargs...) = throw(ArgumentError("Callback parameters not implemented for this spacetime and radiative model"))