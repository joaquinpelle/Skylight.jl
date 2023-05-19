function get_cb_params(model::BosonStarAccretionDisk, configurations)
    
    inner_radius = model.inner_radius
    outer_radius = model.outer_radius

    rmax = get_rmax(configurations)
    
    return BosonStarAccretionDiskCallbackParameters(inner_radius=inner_radius, outer_radius=outer_radius, rmax=rmax)

end

@with_kw struct BosonStarAccretionDiskCallbackParameters <: AbstractCallbackParameters
    
    rmax::Float64
    inner_radius::Float64
    outer_radius::Float64

end

get_callback(::BosonStarAccretionDisk, ::SphericalTopology) = boson_star_accretion_disk_spherical_coordinates_callback()

boson_star_accretion_disk_spherical_coordinates_callback() = VectorContinuousCallback(boson_star_accretion_disk_spherical_coordinates_condition, boson_star_accretion_disk_spherical_coordinates_affect!, 2)

function boson_star_accretion_disk_spherical_coordinates_condition(out, u, t, integrator)
    
    rmax = integrator.p.cb_params.rmax

    out[1] = u[3]-π/2
    out[2] = (rmax-u[2])

end

function boson_star_accretion_disk_spherical_coordinates_affect!(integrator, idx)

    inner_radius = integrator.p.cb_params.inner_radius
    outer_radius = integrator.p.cb_params.outer_radius

    if idx==1 && (inner_radius <= integrator.u[2] <= outer_radius)
        terminate!(integrator)
    elseif idx==2 
        terminate!(integrator) 
    end

end