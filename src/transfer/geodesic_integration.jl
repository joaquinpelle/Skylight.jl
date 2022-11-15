export integrate_geodesics

function integrate_geodesics(initial_data, configurations, cb, cb_params, solver_options)

    N = size(initial_data, 2)  
        
    ensembleprob = set_ensemble_problem(initial_data, configurations, cb_params, solver_options)
    output_data = solve_equations(ensembleprob, cb, solver_options, N)

    return output_data

end

function set_ensemble_problem(initial_data, configurations, cb_params, solver_options)
        
    spacetime = configurations.spacetime
    output_type = solver_options.output_type

    u0 = copy(initial_data[:,1])
    tspan = (0.0, 1e4*cb_params.rmax)
    p = allocate_geodesics_cache(spacetime, cb_params)
    prob = ODEProblem(geodesic_equations!, u0, tspan, p)
    prob_func(prob, i, repeat) = remake(prob, u0 = initial_data[:,i])
    
    return EnsembleProblem(prob; output_func = output_func(output_type), prob_func = prob_func)

end

function solve_equations(ensembleprob, cb, solver_options, N)
    
    method = solver_options.method
    reltol = solver_options.reltol
    abstol = solver_options.abstol
        
    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU

    stats = @timed sim = solve(ensembleprob, method, reltol=reltol, abstol=abstol, callback = cb, EnsembleThreads(); trajectories = N)

    print_stats(stats)

    output_data = collect_output(sim, solver_options.output_type)

    return output_data

end

function collect_output(sim, output_type::SaveEndstate)
    
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


