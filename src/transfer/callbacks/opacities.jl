opacities_callback() = DiscreteCallback(opacities_condition, opacities_affect!)

function opacities_condition(u, t, integrator)
    
    τmax = integrator.p.cb_params.τmax
    NE = integrator.p.NE

    @inbounds for i in 1:NE
        if u[8+i] < τmax return false end
    end

    return true
    
end

opacities_affect!(integrator) = terminate!(integrator)