#Star

function star_cartesian_condition(u,t,integrator)
    
    rmax = integrator.p.geodesic_config.rmax
    rmin = integrator.p.initial_data_config.emission_model.star_radius

    r2 = u[2]^2 + u[3]^2 + u[4]^2
    return (rmax*rmax - r2)*(r2 - rmin*rmin)

end

function star_spherical_condition(u,t,integrator)
    
    rmax = integrator.p.geodesic_config.rmax
    rmin = integrator.p.initial_data_config.emission_model.star_radius

    return (rmax - u[2])*(u[2] - rmin)

end

star_affect!(integrator) = terminate!(integrator)

star_cartesian_callback() = ContinuousCallback(star_cartesian_condition, star_affect!)
star_spherical_callback() = ContinuousCallback(star_spherical_condition, star_affect!)

#Disk cartesian coordinates:

function disk_cartesian_condition(out, u, t, integrator)
    
    rmax = integrator.p.geodesic_config.rmax
    rmin_bound = integrator.p.geodesic_config.rmin_bound
    rmin = integrator.p.initial_data_config.spacetime.rmin+rmin_bound
    
    rho2_A2 = u[2]*u[2] + u[3]*u[3] + u[4]*u[4]-A2
    r2 = 0.5*(rho2_A2+sqrt(rho2_A2*rho2_A2+4.0*A2*u[4]*u[4]))

    out[1] = u[4]
    out[2] = (rmax*rmax-r2)*(r2-rmin*rmin)
    
end

function disk_cartesian_affect!(integrator, idx)
    
    inner_radius = integrator.initial_data_config.emission_model.inner_radius
    outer_radius = integrator.initial_data_config.emission_model.outer_radius

    if idx==1

        rho2_A2 = integrator.u[2]*integrator.u[2] + integrator.u[3]*integrator.u[3] + integrator.u[4]*integrator.u[4]-A2
        r2 = 0.5*(rho2_A2+sqrt(rho2_A2*rho2_A2+4.0*A2*integrator.u[4]*integrator.u[4]))
        if inner_radius*inner_radius <= r2 <= outer_radius*outer_radius terminate!(integrator) end
    
    elseif idx==2 
    
        terminate!(integrator) 
    
    end

end

disk_cartesian_callback() = VectorContinuousCallback(disk_cartesian_condition, disk_cartesian_affect!, 2)

#Disk spherical coordinates

function disk_spherical_condition(out, u, t, integrator)
    
    rmax = integrator.p.geodesic_config.rmax
    rmin_bound = integrator.p.geodesic_config.rmin_bound
    rmin = integrator.p.initial_data_config.spacetime.rmin+rmin_bound

    out[1] = u[3]-π/2
    out[2] = (rmax-u[2])*(u[2]-rmin)

end


function disk_spherical_affect!(integrator, idx)

    inner_radius = integrator.initial_data_config.emission_model.inner_radius
    outer_radius = integrator.initial_data_config.emission_model.outer_radius

    if idx==1 && (inner_radius <= r <= outer_radius)
        terminate!(integrator)
    elseif idx==2 
        terminate!(integrator) 
    end

end

disk_spherical_callback() = VectorContinuousCallback(disk_spherical_condition, disk_spherical_affect!, 2)

#Star across wormhole regular

function star_across_wormhole_condition(out, u, t, integrator)
    
    t,l,θ,φ = u[1:4]

    r = sqrt(l^2+b0^2-Q^2)
    x,y,z = cartesian_from_spherical([r,θ,φ])
    
    l_center =  integrator.initial_data_config.emission_model.l_center
    star_radius =  integrator.initial_data_config.emission_model.star_radius
    
    r_center = sqrt(l_center^2+b0^2-Q^2)
    
    out[1] = (x-r_center)^2 + y^2 + z^2 - star_radius^2  
    out[2] = (Rmax-r)*(Rmax+r)

end

function star_across_wormhole_affect!(integrator,idx)
 
    if idx==1 && integrator.u[2]<0.0  
        terminate!(integrator)
    elseif idx==2 
        terminate!(integrator) 
    end

end

star_across_wormhole_callback() = VectorContinuousCallback(star_across_wormhole_condition, star_across_wormhole_affect!, 2)
