export integrate_transfer

function integrate_transfer(initial_data, configurations, cb, cb_params, solver_options)

    N = size(initial_data, 2)  

    ensembleprob = set_ensemble_transfer_problem(initial_data, configurations, cb_params, solver_options)
    output_data = solve_equations(ensembleprob, cb, solver_options, N)

    return output_data

end

function set_ensemble_transfer_problem(initial_data, configurations, cb_params, solver_options)
        
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    observed_energies = configurations.observed_energies

    output_type = solver_options.output_type

    u0 = copy(initial_data[:,1])
    tspan = (0.0, 1e4*cb_params.rmax)
    p = allocate_transfer_cache(spacetime, model, cb_params, observed_energies)
    prob = ODEProblem(transfer_equations!, u0, tspan, p)
    prob_func(prob, i, repeat) = remake(prob, u0 = initial_data[:,i])
    
    return EnsembleProblem(prob; output_func = output_func(output_type), prob_func = prob_func)

end