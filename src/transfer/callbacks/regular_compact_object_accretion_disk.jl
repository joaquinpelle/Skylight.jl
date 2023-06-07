function callback_parameters(::AbstractRegularCompactObjectSpacetime, model::AbstractAccretionDisk, configurations)
    inner_radius = model.inner_radius
    outer_radius = model.outer_radius
    rmax = max_radius(configurations)
    return RegularCompactObjectAccretionDiskCallbackParameters(inner_radius=inner_radius, outer_radius=outer_radius, rmax=rmax)
end

@with_kw struct RegularCompactObjectAccretionDiskCallbackParameters <: AbstractCallbackParameters
    rmax::Float64
    inner_radius::Float64
    outer_radius::Float64

    @assert inner_radius >= 0.0 "Inner radius must be non-negative"
    @assert outer_radius >= inner_radius "Outer radius must be larger than inner radius"
    @assert rmax > outer_radius "rmax must be larger than outer_radius"
end

callback(::AbstractRegularCompactObjectSpacetime, ::AbstractAccretionDisk, ::SphericalTopology) = regular_spacetime_accretion_disk_spherical_coordinates_callback()

regular_spacetime_accretion_disk_spherical_coordinates_callback() = VectorContinuousCallback(regular_spacetime_accretion_disk_spherical_coordinates_condition, regular_spacetime_accretion_disk_spherical_coordinates_affect!, 2)

function regular_spacetime_accretion_disk_spherical_coordinates_condition(out, u, t, integrator)    
    rmax = integrator.p.cbp.rmax
    out[1] = u[3]-π/2
    out[2] = (rmax-u[2])
end

function regular_spacetime_accretion_disk_spherical_coordinates_affect!(integrator, idx)
    inner_radius = integrator.p.cbp.inner_radius
    outer_radius = integrator.p.cbp.outer_radius

    if (inner_radius <= integrator.u[2] <= outer_radius) || idx == 2
        terminate!(integrator) 
    end
end