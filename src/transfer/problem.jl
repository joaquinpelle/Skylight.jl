integrate(initial_data, configurations::AbstractConfigurations, cb, cbp; kwargs...) = integrate(isvacuum(configurations), initial_data, configurations, cb, cbp; kwargs...)
ensemble_problem(initial_data::AbstractMatrix, configurations::AbstractConfigurations, cbp::AbstractCallbackParameters, args...) = ensemble_problem(isvacuum(configurations), initial_data, configurations, cbp, args...)

function integrate(::Vacuum, initial_data, configurations, cb, cbp; method = VCABM(), kwargs...)
    N = size(initial_data, 2)  
    ensembleprob = ensemble_problem(initial_data, configurations, cbp)
    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU
    stats = @timed sim = solve(ensembleprob, method, EnsembleThreads(); callback = cb, trajectories = N, kwargs...)
    print_stats(stats)
    return collect_run(sim, cb, cbp, method; kwargs...)
end

function integrate(::NonVacuum, initial_data, configurations, cb, cbp; τmax=2.0, method=VCABM(), kwargs...)
    N = size(initial_data, 2)  
    ensembleprob = ensemble_problem(initial_data, configurations, cbp, τmax)
    full_cb = CallbackSet(cb, opacities_callback())  
    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU
    stats = @timed sim = solve(ensembleprob, method, EnsembleThreads(); callback = full_cb, trajectories = N, kwargs...)
    print_stats(stats)
    return collect_run(sim, cb, cbp, τmax, method; kwargs...)
end

function ensemble_problem(::Vacuum, initial_data::AbstractMatrix, configurations::AbstractConfigurations, cbp::AbstractCallbackParameters)
    u0 = SVector{8, Float64}(initial_data[:,1]...)
    tspan = (0.0, 1e4*cbp.rmax)
    p = transfer_cache(configurations, cbp)
    prob = ODEProblem(equations(configurations), u0, tspan, p)
    output_func(sol, i) = (sol[end], false)
    prob_func(prob, i, repeat) = remake(prob, u0 = SVector{8,Float64}(initial_data[:,i]...))
    return EnsembleProblem(prob; output_func = output_func, prob_func = prob_func)
end

function ensemble_problem(::NonVacuum, initial_data::AbstractMatrix, configurations::AbstractConfigurations, cbp::AbstractCallbackParameters, τmax::Real)
    NE = length(configurations.observation_energies)
    Nvars = 8+2*NE
    u0 = SVector{Nvars, Float64}(initial_data[:,1]...)
    tspan = (0.0, 1e4*cbp.rmax)
    p = transfer_cache(configurations, cbp, τmax)
    prob = ODEProblem(equations(configurations), u0, tspan, p)
    output_func(sol, i) = (sol[end], false)
    prob_func(prob, i, repeat) = remake(prob, u0 = SVector{Nvars,Float64}(initial_data[:,i]...))
    return EnsembleProblem(prob; output_func = output_func, prob_func = prob_func)
end

function collect_run(sim, cb, cbp, args...; kwargs...)
    output_data = collect_output(sim)
    kwargs_dict = collect_args(args...; kwargs...)
    return Run(output_data, cb, cbp, kwargs_dict)
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

function collect_args(args...; kwargs...)
    args_dict = Dict{Symbol, Any}([Symbol("arg", i) => arg for (i, arg) in enumerate(args)]...)
    kwargs_dict = Dict{Symbol, Any}(pairs(kwargs))
    merged_dict = merge(args_dict, kwargs_dict)
    return merged_dict
end

function print_stats(stats)
    println("Equations integration stats:")
    println("Wall clock time: ",stats.time," seconds")
    println("Memory allocated: ",Base.format_bytes(stats.bytes))
end

to_tuple(run::Run) = (run.output_data, run.callback, run.callback_parameters, run.kwargs)