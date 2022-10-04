function dump_vector_in(container, vector)
    container[:,5] = vector
end

function dump_∂t_in!(container)
    container[:,5] = ∂t()
end

function dump_metric_in!(container,position,spacetime::Spacetime)

    pars = spacetime.parameters
    metric! = spacetime.metric!
    
    @views g = container[:,1:4]
    metric!(g,position,pars)
    
end