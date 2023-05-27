function get_cb_params(::AbstractSpacetime, model::SyntheticPolarCap, configurations) 
    rmin = model.star_radius
    rmax = get_rmax(configurations)
    return NeutronStarHotSpotsCallbackParameters(rmin=rmin, rmax=rmax)
end

@with_kw struct NeutronStarHotSpotsCallbackParameters <: AbstractCallbackParameters
    rmax::Float64
    rmin::Float64
end

get_callback(::AbstractSpacetime, ::SyntheticPolarCap, ::CartesianTopology) = star_cartesian_coordinates_callback()
get_callback(::AbstractSpacetime, ::SyntheticPolarCap, ::SphericalTopology) = star_spherical_coordinates_callback()

star_cartesian_coordinates_callback() = ContinuousCallback(star_cartesian_coordinates_condition, star_affect!)
star_spherical_coordinates_callback() = ContinuousCallback(star_spherical_coordinates_condition, star_affect!)

function star_cartesian_coordinates_condition(u,t,integrator)
    
    rmax = integrator.p.cb_params.rmax
    rmin = integrator.p.cb_params.rmin

    r2 = u[2]^2 + u[3]^2 + u[4]^2
    return (rmax*rmax - r2)*(r2 - rmin*rmin)

end

function star_spherical_coordinates_condition(u,t,integrator)
    
    rmax = integrator.p.cb_params.rmax
    rmin = integrator.p.cb_params.rmin

    return (rmax - u[2])*(u[2] - rmin)

end

star_affect!(integrator) = terminate!(integrator)

