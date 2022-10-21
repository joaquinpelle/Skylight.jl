@with_kw mutable struct ThreadCache{T<:ChristoffelCache}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T

end

@with_kw mutable struct GeodesicsCache{S<:Spacetime,M<:EmissionModel,C<:CoordinateSystemKind,T<:ChristoffelCache}
    spacetime::S
    model::M
    multi_thread::Array{ThreadCache{T},1}
    coord_system::C = coordinate_system_kind(spacetime)
end

function allocate_single_thread_cache(spacetime) 
    return ThreadCache(christoffel_cache = allocate_christoffel_cache(spacetime))
end

function allocate_multi_thread_cache(spacetime)
    return [allocate_single_thread_cache(spacetime) for i in 1:Threads.nthreads()]
end

function allocate_geodesics_cache(spacetime)
    return GeodesicsCache(spacetime=spacetime, multi_thread=allocate_multi_thread_cache(spacetime))
end