function transfer_cache(configurations::AbstractConfigurations, args...)
    transfer_cache(isvacuum(configurations), configurations, args...)
end

# Vacuum
function transfer_cache(::Vacuum, configurations, cbp)
    spacetime = configurations.spacetime
    return VacuumCache(spacetime, cbp, vacuum_multi_thread_cache(spacetime))
end

function vacuum_multi_thread_cache(spacetime)
    return [vacuum_single_thread_cache(spacetime) for i in 1:maxthreadid()]
end

function vacuum_single_thread_cache(spacetime)
    return VacuumThreadCache(spacetime_cache = allocate_cache(spacetime),
        christoffel_cache = allocate_christoffel_cache(spacetime))
end

# Non vacuum
function transfer_cache(::NonVacuum, configurations, cbp, τmax)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    observation_energies = configurations.observation_energies
    NE = length(observation_energies)
    return NonVacuumCache(spacetime,
        model,
        coordinates_topology(spacetime),
        cbp,
        τmax,
        observation_energies,
        geometrized_to_CGS(1.0, Dimensions.length, configurations),
        length(observation_energies),
        non_vacuum_multi_thread_cache(spacetime, model, NE))
end

function non_vacuum_multi_thread_cache(spacetime, model, NE)
    return [non_vacuum_single_thread_cache(spacetime, model, NE) for i in 1:maxthreadid()]
end

function non_vacuum_single_thread_cache(spacetime, model, NE)
    return NonVacuumThreadCache(spacetime_cache = allocate_cache(spacetime),
        model_cache = allocate_cache(model),
        christoffel_cache = allocate_christoffel_cache(spacetime),
        ε = zeros(NE),
        αε = zeros(NE),
        jε = zeros(NE))
end
