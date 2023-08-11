function callback_parameters(spacetime::AbstractBlackHoleSpacetime, ::IonTorus, configurations; rhorizon_bound)
    rhorizon = event_horizon_radius(spacetime)
    rmin = rhorizon + rhorizon_bound
    rmax = max_radius(configurations)
    return IonTorusCallbackParameters(rmin=rmin, rmax=rmax)
end

@with_kw struct IonTorusCallbackParameters <: AbstractCallbackParameters
    rmin::Float64
    rmax::Float64
    @assert 0.0 ≤ rmin < rmax "rmin must be smaller than rmax and non-negative"
end

callback(::AbstractBlackHoleSpacetime, ::IonTorus, ::SphericalTopology) = ContinuousCallback(iontorus_condition, iontorus_affect!)

function iontorus_condition(u, t, integrator)
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    return (rmax-u[2])*(u[2]-rmin)
end

iontorus_affect!(integrator) = terminate!(integrator)