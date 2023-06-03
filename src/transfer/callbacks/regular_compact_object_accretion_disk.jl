function get_cb_params(::AbstractRegularCompactObjectSpacetime, model::AbstractAccretionDisk, configurations)
    inner_radius = model.inner_radius
    outer_radius = model.outer_radius
    rmax = get_rmax(configurations)
    return RegularCompactObjectAccretionDiskCallbackParameters(inner_radius=inner_radius, outer_radius=outer_radius, rmax=rmax)
end

@with_kw struct RegularCompactObjectAccretionDiskCallbackParameters <: AbstractCallbackParameters
    rmax::Float64
    inner_radius::Float64
    outer_radius::Float64
end

get_callback(::AbstractRegularCompactObjectSpacetime, ::AbstractAccretionDisk, ::SphericalTopology) = regular_spacetime_accretion_disk_spherical_coordinates_callback()

regular_spacetime_accretion_disk_spherical_coordinates_callback() = VectorContinuousCallback(regular_spacetime_accretion_disk_spherical_coordinates_condition, regular_spacetime_accretion_disk_spherical_coordinates_affect!, 2)

function regular_spacetime_accretion_disk_spherical_coordinates_condition(out, u, t, integrator)    
    rmax = integrator.p.cb_params.rmax
    out[1] = u[3]-π/2
    out[2] = (rmax-u[2])
end

function regular_spacetime_accretion_disk_spherical_coordinates_affect!(integrator, idx)
    inner_radius = integrator.p.cb_params.inner_radius
    outer_radius = integrator.p.cb_params.outer_radius

    if (inner_radius <= integrator.u[2] <= outer_radius) || idx == 2
        terminate!(integrator) 
    end
end