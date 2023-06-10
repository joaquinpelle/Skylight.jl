function callback_parameters(spacetime::AbstractBlackHoleSpacetime, model::AbstractAccretionDisk, configurations; rhorizon_bound)
    inner_radius = model.inner_radius
    outer_radius = model.outer_radius

    rhorizon = event_horizon_radius(spacetime)
    rmin = rhorizon + rhorizon_bound
    rmax = max_radius(configurations)
    
    return BlackHoleAccretionDiskCallbackParameters(inner_radius=inner_radius, outer_radius=outer_radius, rmin=rmin, rmax=rmax)
end

@with_kw struct BlackHoleAccretionDiskCallbackParameters <: AbstractCallbackParameters
    rmin::Float64
    rmax::Float64
    inner_radius::Float64
    outer_radius::Float64
    
    @assert inner_radius >= 0.0 "Inner radius must be non-negative"
    @assert outer_radius >= inner_radius "Outer radius must be larger than inner radius"
    @assert rmin <= inner_radius "rmin must be smaller than inner radius"
    @assert rmax > outer_radius "rmax must be larger than outer_radius"
end

callback(::AbstractBlackHoleSpacetime, ::AbstractAccretionDisk, ::SphericalTopology) = black_hole_accretion_disk_spherical_coordinates_callback()
callback(::SchwarzschildSpacetimeKerrSchildCoordinates, ::AbstractAccretionDisk, ::CartesianTopology) = nonspinning_black_hole_accretion_disk_cartesian_coordinates_callback()
callback(::KerrSpacetimeKerrSchildCoordinates, ::AbstractAccretionDisk, ::CartesianTopology) = spinning_black_hole_accretion_disk_cartesian_coordinates_callback()

spinning_black_hole_accretion_disk_cartesian_coordinates_callback() = VectorContinuousCallback(spinning_black_hole_accretion_disk_cartesian_coordinates_condition, spinning_black_hole_accretion_disk_cartesian_coordinates_affect!, 2)
nonspinning_black_hole_accretion_disk_cartesian_coordinates_callback() = VectorContinuousCallback(nonspinning_black_hole_accretion_disk_cartesian_coordinates_condition, nonspinning_black_hole_accretion_disk_cartesian_coordinates_affect!, 2)
black_hole_accretion_disk_spherical_coordinates_callback() = VectorContinuousCallback(black_hole_accretion_disk_spherical_coordinates_condition, black_hole_accretion_disk_spherical_coordinates_affect!, 2)

function nonspinning_black_hole_accretion_disk_cartesian_coordinates_condition(out, u, t, integrator)
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    r2 = u[2]*u[2] + u[3]*u[3] + u[4]*u[4]
    out[1] = u[4]
    out[2] = (rmax*rmax-r2)*(r2-rmin*rmin)
end

function nonspinning_black_hole_accretion_disk_cartesian_coordinates_affect!(integrator, idx)
    inner_radius = integrator.p.cbp.inner_radius
    outer_radius = integrator.p.cbp.outer_radius
    if idx==1
        r2 = integrator.u[2]*integrator.u[2] + integrator.u[3]*integrator.u[3] + integrator.u[4]*integrator.u[4]
        if inner_radius*inner_radius <= r2 <= outer_radius*outer_radius terminate!(integrator) end
    elseif idx==2 
        terminate!(integrator)
    end
end

function spinning_black_hole_accretion_disk_cartesian_coordinates_condition(out, u, t, integrator)
    a = integrator.p.spacetime.a
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    rho2_a2 = u[2]*u[2] + u[3]*u[3] + u[4]*u[4]-a*a
    r2 = 0.5*(rho2_a2+sqrt(rho2_a2*rho2_a2+4.0*a*a*u[4]*u[4]))
    out[1] = u[4]
    out[2] = (rmax*rmax-r2)*(r2-rmin*rmin)
end

function spinning_black_hole_accretion_disk_cartesian_coordinates_affect!(integrator, idx)
    a = integrator.p.spacetime.a
    inner_radius = integrator.p.cbp.inner_radius
    outer_radius = integrator.p.cbp.outer_radius
    if idx==1
        rho2_a2 = integrator.u[2]*integrator.u[2] + integrator.u[3]*integrator.u[3] + integrator.u[4]*integrator.u[4]-a*a
        r2 = 0.5*(rho2_a2+sqrt(rho2_a2*rho2_a2+4.0*a*a*integrator.u[4]*integrator.u[4]))
        if inner_radius*inner_radius <= r2 <= outer_radius*outer_radius terminate!(integrator) end
    elseif idx==2 
        terminate!(integrator)
    end
end

function black_hole_accretion_disk_spherical_coordinates_condition(out, u, t, integrator)
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    out[1] = u[3]-π/2
    out[2] = (rmax-u[2])*(u[2]-rmin)
end

function black_hole_accretion_disk_spherical_coordinates_affect!(integrator, idx)
    inner_radius = integrator.p.cbp.inner_radius
    outer_radius = integrator.p.cbp.outer_radius
    if (inner_radius <= integrator.u[2] <= outer_radius) || idx == 2
        terminate!(integrator) 
    end
end