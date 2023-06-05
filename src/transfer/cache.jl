# Vacuum

function allocate_cache(configurations::VacuumConfigurations, cbp)
    spacetime = configurations.spacetime
    return VacuumCache(spacetime, cbp, allocate_vacuum_multi_thread_cache(spacetime))
end

function allocate_vacuum_multi_thread_cache(spacetime)
    return [allocate_vacuum_single_thread_cache(spacetime) for i in 1:Threads.nthreads()]
end

function allocate_vacuum_single_thread_cache(spacetime) 
    return VacuumThreadCache(christoffel_cache = allocate_christoffel_cache(spacetime))
end

# Non vacuum

function allocate_cache(configurations::NonVacuumConfigurations, cbp, τmax)

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    observation_energies = configurations.observation_energies
    NE = length(observation_energies)

    return NonVacuumCache(spacetime, 
                         model, 
                         cbp,
                         τmax, 
                         observation_energies,
                         length(observation_energies), 
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