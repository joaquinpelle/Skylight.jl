function callback_parameters(::AbstractSpacetime, model::SyntheticPolarCap, configurations) 
    rmin = model.star_radius
    rmax = max_radius(configurations)
    return NeutronStarHotSpotsCallbackParameters(rmin=rmin, rmax=rmax)
end

@with_kw struct NeutronStarHotSpotsCallbackParameters <: AbstractCallbackParameters
    rmax::Float64
    rmin::Float64

    @assert rmin >= 0.0 "rmin must be non-negative"
    @assert rmax > rmin "rmax must be larger than rmin"
end

callback(::AbstractSpacetime, ::SyntheticPolarCap, ::CartesianTopology) = star_cartesian_coordinates_callback()
callback(::AbstractSpacetime, ::SyntheticPolarCap, ::SphericalTopology) = star_spherical_coordinates_callback()

star_cartesian_coordinates_callback() = ContinuousCallback(star_cartesian_coordinates_condition, terminate!)
star_spherical_coordinates_callback() = ContinuousCallback(star_spherical_coordinates_condition, terminate!)

function star_cartesian_coordinates_condition(u,t,integrator)
    
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin

    r2 = u[2]^2 + u[3]^2 + u[4]^2
    return (rmax*rmax - r2)*(r2 - rmin*rmin)
end

function star_spherical_coordinates_condition(u,t,integrator)
    
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin

    return (rmax - u[2])*(u[2] - rmin)
end

