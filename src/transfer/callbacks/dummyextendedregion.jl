# Dummy extended region

function callback_parameters(::AbstractBlackHoleSpacetime, ::DummyExtendedRegion, configurations; rhorizon_bound) 
    rmax = max_radius(configurations)
    rhorizon = event_horizon_radius(configurations.spacetime)
    rmin = rhorizon + rhorizon_bound
    return DummyExtendedRegionCallbackParameters(rmin=rmin, rmax=rmax)
end

@with_kw struct DummyExtendedRegionCallbackParameters <: AbstractCallbackParameters    
    rmin::Float64
    rmax::Float64
end

callback(::AbstractSpacetime, ::DummyExtendedRegion, ::AbstractCoordinatesTopology) = dummy_extended_region_callback()

dummy_extended_region_callback() = ContinuousCallback(dummy_extended_region_condition, dummy_extended_region_affect!, 2)

function dummy_extended_region_condition(u, t, integrator)    
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    return (rmax*-u[2])*(u[2]-rmin)
end

dummy_extended_region_affect!(integrator) = terminate!(integrator)