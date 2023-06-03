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