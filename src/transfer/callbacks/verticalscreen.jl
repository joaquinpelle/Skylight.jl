@with_kw struct VerticalScreenCallbackParameters <: AbstractCallbackParameters
    rmin::Float64
    rmax::Float64
    x::Float64
    @assert rmax > 0.0 "rmax must be positive"
end

function callback_parameters(spacetime::AbstractBlackHoleSpacetime, model::VerticalScreen, configurations; rmax, rhorizon_bound) 
    x = model.x
    rhorizon = event_horizon_radius(spacetime)
    rmin = rhorizon + rhorizon_bound
    return VerticalScreenCallbackParameters(rmin=rmin, rmax=rmax, x=x)
end

callback(::AbstractBlackHoleSpacetime, ::VerticalScreen, ::CartesianTopology) = screen_behind_black_hole_callback()

screen_behind_black_hole_callback() = VectorContinuousCallback(screen_behind_black_hole_condition, screen_behind_black_hole_affect!, 2)

function screen_behind_black_hole_condition(out, u, t, integrator)
    a = integrator.p.spacetime.a
    xscreen = integrator.p.cbp.x
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    rho2_a2 = u[2]*u[2] + u[3]*u[3] + u[4]*u[4]-a*a
    r2 = 0.5*(rho2_a2+sqrt(rho2_a2*rho2_a2+4.0*a*a*u[4]*u[4]))
    out[1] = u[2]-xscreen
    out[2] = (rmax*rmax-r2)*(r2-rmin*rmin)
end

screen_behind_black_hole_affect!(integrator,idx) = terminate!(integrator)