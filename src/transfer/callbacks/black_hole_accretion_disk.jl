function get_cb_params(spacetime::AbstractBlackHoleSpacetime, model::AbstractAccretionDisk, configurations; rhorizon_bound)
    inner_radius = model.inner_radius
    outer_radius = model.outer_radius

    rhorizon = event_horizon_radius(spacetime)
    rmin = rhorizon + rhorizon_bound
    rmax = get_rmax(configurations)
    
    return BlackHoleAccretionDiskCallbackParameters(inner_radius=inner_radius, outer_radius=outer_radius, rmin=rmin, rmax=rmax)
end

@with_kw struct BlackHoleAccretionDiskCallbackParameters <: AbstractCallbackParameters
    rmax::Float64
    rmin::Float64
    inner_radius::Float64
    outer_radius::Float64
end

get_callback(::AbstractBlackHoleSpacetime, ::AbstractAccretionDisk, ::CartesianTopology) = black_hole_accretion_disk_cartesian_coordinates_callback()
get_callback(::AbstractBlackHoleSpacetime, ::AbstractAccretionDisk, ::SphericalTopology) = black_hole_accretion_disk_spherical_coordinates_callback()

black_hole_accretion_disk_cartesian_coordinates_callback() = VectorContinuousCallback(black_hole_accretion_disk_cartesian_coordinates_condition, black_hole_accretion_disk_cartesian_coordinates_affect!, 2)
black_hole_accretion_disk_spherical_coordinates_callback() = VectorContinuousCallback(black_hole_accretion_disk_spherical_coordinates_condition, black_hole_accretion_disk_spherical_coordinates_affect!, 2)

function black_hole_accretion_disk_cartesian_coordinates_condition(out, u, t, integrator)
    
    a = integrator.p.spacetime.a

    rmax = integrator.p.cb_params.rmax
    rmin = integrator.p.cb_params.rmin
    
    rho2_a2 = u[2]*u[2] + u[3]*u[3] + u[4]*u[4]-a*a
    r2 = 0.5*(rho2_a2+sqrt(rho2_a2*rho2_a2+4.0*a*a*u[4]*u[4]))

    out[1] = u[4]
    out[2] = (rmax*rmax-r2)*(r2-rmin*rmin)
end

function black_hole_accretion_disk_cartesian_coordinates_affect!(integrator, idx)
    
    a = integrator.p.spacetime.a
    inner_radius = integrator.p.cb_params.inner_radius
    outer_radius = integrator.p.cb_params.outer_radius

    if idx==1
        rho2_a2 = integrator.u[2]*integrator.u[2] + integrator.u[3]*integrator.u[3] + integrator.u[4]*integrator.u[4]-a*a
        r2 = 0.5*(rho2_a2+sqrt(rho2_a2*rho2_a2+4.0*a*a*integrator.u[4]*integrator.u[4]))
        if inner_radius*inner_radius <= r2 <= outer_radius*outer_radius terminate!(integrator) end
    elseif idx==2 
        terminate!(integrator)
    end
end

function black_hole_accretion_disk_spherical_coordinates_condition(out, u, t, integrator)
    
    rmax = integrator.p.cb_params.rmax
    rmin = integrator.p.cb_params.rmin

    out[1] = u[3]-π/2
    out[2] = (rmax-u[2])*(u[2]-rmin)

end

function black_hole_accretion_disk_spherical_coordinates_affect!(integrator, idx)

    inner_radius = integrator.p.cb_params.inner_radius
    outer_radius = integrator.p.cb_params.outer_radius

    if idx==1 && (inner_radius <= integrator.u[2] <= outer_radius)
        terminate!(integrator)
    elseif idx==2 
        terminate!(integrator) 
    end
end