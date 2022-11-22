# Dummy extended region

function get_cb_params(model::DummyExtendedRegion, configurations; rhorizon_bound) 

    rmax = get_rmax(configurations)

    rhorizon = event_horizon_radius(configurations.spacetime)
    rmin = rhorizon + rhorizon_bound

    return DummyExtendedRegionCallbackParameters(rmin=rmin, rmax=rmax)

end

@with_kw struct DummyExtendedRegionCallbackParameters <: CallbackParameters
    
    rmin::Float64
    rmax::Float64

end

get_callback(model::DummyExtendedRegion, coord_system) = dummy_extended_region_callback()

dummy_extended_region_callback() = ContinuousCallback(dummy_extended_region_condition, dummy_extended_region_affect!, 2)

function dummy_extended_region_condition(u, t, integrator)
    
    rmax = integrator.p.cb_params.rmax
    rmin = integrator.p.cb_params.rmin

    return (rmax*-u[2])*(u[2]-rmin)
    
end

dummy_extended_region_affect!(integrator) = terminate!(integrator)