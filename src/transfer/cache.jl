# Vacuum

function allocate_cache(configurations::VacuumConfigurations, cb_params)
    spacetime = configurations.spacetime
    return VacuumCache(spacetime, cb_params, allocate_vacuum_multi_thread_cache(spacetime))
end

function allocate_vacuum_multi_thread_cache(spacetime)
    return [allocate_vacuum_single_thread_cache(spacetime) for i in 1:Threads.nthreads()]
end

function allocate_vacuum_single_thread_cache(spacetime) 
    return VacuumThreadCache(christoffel_cache = allocate_christoffel_cache(spacetime))
end


@with_kw mutable struct VacuumThreadCache{T}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T

end

mutable struct VacuumCache{S, C, T}
    
    spacetime::S
    cb_params::C
    multi_thread::Array{VacuumThreadCache{T},1}

end

# Non vacuum

function allocate_cache(configurations::NonVacuumConfigurations, cb_params, τmax)

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    observed_energies = configurations.observed_energies
    NE = length(observed_energies)

    return NonVacuumCache(spacetime, 
                         model, 
                         cb_params,
                         τmax, 
                         observed_energies,
                         length(observed_energies), 
                         allocate_non_vacuum_multi_thread_cache(spacetime, NE))

end

function allocate_non_vacuum_multi_thread_cache(spacetime, NE)
   
    return [allocate_non_vacuum_single_thread_cache(spacetime, NE) for i in 1:Threads.nthreads()]

end

function allocate_non_vacuum_single_thread_cache(spacetime, NE)
    
    return NonVacuumThreadCache(christoffel_cache = allocate_christoffel_cache(spacetime),
                               ε = zeros(NE),
                               αε = zeros(NE),
                               jε = zeros(NE))

end

@with_kw mutable struct NonVacuumThreadCache{T}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T
    ε::Array{Float64, 1}
    αε::Array{Float64, 1}
    jε::Array{Float64, 1}

end

mutable struct NonVacuumCache{S<:Spacetime, M<:RadiativeModel, C<:CallbackParameters, T<:ChristoffelCache}
    
    spacetime::S
    model::M
    cb_params::C
    τmax::Float64
    observed_energies::Array{Float64, 1}
    NE::Int64
    multi_thread::Array{NonVacuumThreadCache{T},1}

end