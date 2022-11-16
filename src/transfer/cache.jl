@with_kw mutable struct GeodesicThreadCache{T<:ChristoffelCache}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T

end

mutable struct GeodesicsCache{S<:Spacetime, C<:CallbackParameters, T<:ChristoffelCache}
    
    spacetime::S
    cb_params::C
    multi_thread::Array{GeodesicThreadCache{T},1}

end

function allocate_geodesics_cache(spacetime, cb_params)
    return GeodesicsCache(spacetime, cb_params, allocate_geodesics_multi_thread_cache(spacetime))
end

function allocate_geodesics_multi_thread_cache(spacetime)
    return [allocate_geodesic_single_thread_cache(spacetime) for i in 1:Threads.nthreads()]
end

function allocate_geodesic_single_thread_cache(spacetime) 
    return GeodesicThreadCache(christoffel_cache = allocate_christoffel_cache(spacetime))
end


@with_kw mutable struct TransferThreadCache{T<:ChristoffelCache}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T
    ε::Array{Float64, 1}
    αε::Array{Float64, 1}
    jε::Array{Float64, 1}

end

mutable struct TransferCache{S<:Spacetime, M<:RadiativeModel, C<:CallbackParameters, T<:ChristoffelCache}
    
    spacetime::S
    model::M
    cb_params::C
    observed_energies::Array{Float64, 1}
    NE::Int64
    multi_thread::Array{TransferThreadCache{T},1}

end

function allocate_transfer_cache(spacetime, model, cb_params, observed_energies)

    NE = length(observed_energies)

    return TransferCache(spacetime, 
                         model, 
                         cb_params, 
                         observed_energies,
                         NE, 
                         allocate_transfer_multi_thread_cache(spacetime, NE))

end

function allocate_transfer_multi_thread_cache(spacetime, NE)
   
    return [allocate_transfer_single_thread_cache(spacetime, NE) for i in 1:Threads.nthreads()]

end

function allocate_transfer_single_thread_cache(spacetime, NE)
    
    return TransferThreadCache(christoffel_cache = allocate_christoffel_cache(spacetime),
                               ε = zeros(NE),
                               αε = zeros(NE),
                               jε = zeros(NE))

end