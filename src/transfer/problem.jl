export integrate

function integrate(initial_data, configurations::VacuumConfigurations, cb, cb_params; method = VCABM(), kwargs...)

    N = size(initial_data, 2)  

    ensembleprob = set_ensemble_problem(initial_data, configurations, cb_params)

    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU
    stats = @timed sim = solve(ensembleprob, method, EnsembleThreads(); callback = cb, trajectories = N, kwargs...)
    
    print_stats(stats)

    output_data = collect_output(sim)

    return output_data

end

function integrate(initial_data, configurations::NonVacuumConfigurations, cb, cb_params; τmax=2.0, method=VCABM(), kwargs...)

    N = size(initial_data, 2)  

    ensembleprob = set_ensemble_problem(initial_data, configurations, cb_params, τmax)

    full_cb = CallbackSet(cb, opacities_callback())  

    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU
    stats = @timed sim = solve(ensembleprob, method, EnsembleThreads(); callback = full_cb, trajectories = N, kwargs...)
    
    print_stats(stats)

    output_data = collect_output(sim)

    return output_data

end

function set_ensemble_problem(initial_data, configurations::VacuumConfigurations, cb_params)
        
    u0 = copy(initial_data[:,1])
    tspan = (0.0, 1e4*cb_params.rmax)
    p = allocate_cache(configurations, cb_params)
    prob = ODEProblem(equations(configurations), u0, tspan, p)

    output_func(sol, i) = (sol[end], false)
    prob_func(prob, i, repeat) = remake(prob, u0 = initial_data[:,i])

    return EnsembleProblem(prob; output_func = output_func, prob_func = prob_func)
    
end

function set_ensemble_problem(initial_data, configurations::NonVacuumConfigurations, cb_params, τmax)
        
    u0 = copy(initial_data[:,1])
    tspan = (0.0, 1e4*cb_params.rmax)
    p = allocate_cache(configurations, cb_params, τmax)
    prob = ODEProblem(equations(configurations), u0, tspan, p)

    output_func(sol, i) = (sol[end], false)
    prob_func(prob, i, repeat) = remake(prob, u0 = initial_data[:,i])

    return EnsembleProblem(prob; output_func = output_func, prob_func = prob_func)
    
end

function collect_output(sim)
    
    N = size(sim.u, 1)
    Nvars = length(sim.u[1])

    output_data = zeros(Nvars,N)
    
    for i in 1:N
        output_data[:,i] = sim.u[i]
    end

    return output_data

end


function print_stats(stats)
   
    println("Solve stats:")
    println("Wall clock time: ",stats.time," seconds")
    println("Memory allocated: ",Base.format_bytes(stats.bytes))

end
