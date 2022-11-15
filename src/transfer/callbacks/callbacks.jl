get_callback(model::NeutronStarHotSpots, coord_system::CartesianClass) = star_cartesian_coordinates_callback()
get_callback(model::NeutronStarHotSpots, coord_system::SphericalClass) = star_spherical_coordinates_callback()
get_callback(model::BlackHoleAccretionDisk, coord_system::CartesianClass) = black_hole_accretion_disk_cartesian_coordinates_callback()
get_callback(model::BlackHoleAccretionDisk, coord_system::SphericalClass) = black_hole_accretion_disk_spherical_coordinates_callback()
get_callback(model::StarAcrossWormhole, coord_system) = star_across_wormhole_callback()
get_callback(model::DummyExtendedRegion, coord_system) = dummy_extended_region_callback()

#Neutron star hot spots

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


#Black hole accretion disk

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

#Star across wormhole

star_across_wormhole_callback() = VectorContinuousCallback(star_across_wormhole_condition, star_across_wormhole_affect!, 2)

function star_across_wormhole_condition(out, u, t, integrator)
    
    l_center = integrator.p.cb_params.l_center
    star_radius = integrator.p.cb_params.star_radius
    rmax = integrator.p.cb_params.rmax

    b0 = integrator.p.spacetime.b0
    Q = integrator.p.spacetime.Q

    t,l,θ,φ = u[1:4]

    r = sqrt(l^2+b0^2-Q^2)

    x,y,z = cartesian_from_spherical([r,θ,φ])
    
    r_center = sqrt(l_center^2+b0^2-Q^2)
    
    out[1] = (x-r_center)^2 + y^2 + z^2 - star_radius^2  
    out[2] = (rmax-r)*(rmax+r)

end

function star_across_wormhole_affect!(integrator,idx)
 
    if idx==1 && integrator.u[2]<0.0  
        terminate!(integrator)
    elseif idx==2 
        terminate!(integrator) 
    end

end

# Dummy extended region

dummy_extended_region_callback() = CallbackSet(ContinuousCallback(dummy_extended_region_condition, dummy_extended_region_affect!, 2),
                                               opacities_callback())


function dummy_extended_region_condition(u, t, integrator)
    
    rmax = integrator.p.cb_params.rmax
    rmin = integrator.p.cb_params.rmin

    return (rmax*-u[2])*(u[2]-rmin)
    
end

dummy_extended_region_affect!(integrator) = terminate!(integrator)

# Opacities callback

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