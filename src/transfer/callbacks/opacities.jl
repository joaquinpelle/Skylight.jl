opacities_callback() = DiscreteCallback(opacities_condition, terminate!)

function opacities_condition(u, t, integrator)
    τmax = integrator.p.τmax
    NE = integrator.p.NE

    @inbounds begin
        for i in 1:NE
            if u[8 + i] < τmax
                return false
            end
        end
    end
    return true
end
