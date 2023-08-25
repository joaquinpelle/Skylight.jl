@with_kw struct StarAcrossWormholeCallbackParameters <: AbstractCallbackParameters
    rmax::Float64
    l_center::Float64
    star_radius::Float64

    @assert rmax>0.0 "rmax must be positive"
    @assert star_radius>0.0 "star_radius must be positive"
end

function callback_parameters(::ChargedWormholeSpacetimeRegularCoordinates,
    model::StarAcrossWormhole,
    configurations)
    l_center = model.l_center
    star_radius = model.star_radius
    rmax = max_radius(configurations)

    return StarAcrossWormholeCallbackParameters(l_center = l_center,
        star_radius = star_radius,
        rmax = rmax)
end

function callback(::ChargedWormholeSpacetimeRegularCoordinates,
    ::StarAcrossWormhole,
    ::AbstractCoordinatesTopology)
    star_across_wormhole_callback()
end

function star_across_wormhole_callback()
    VectorContinuousCallback(star_across_wormhole_condition,
        star_across_wormhole_affect!,
        2)
end

function star_across_wormhole_condition(out, u, t, integrator)
    l_center = integrator.p.cbp.l_center
    star_radius = integrator.p.cbp.star_radius
    rmax = integrator.p.cbp.rmax

    b0 = integrator.p.spacetime.b0
    Q = integrator.p.spacetime.Q

    t, l, θ, φ = u[1:4]

    r = sqrt(l^2 + b0^2 - Q^2)

    x, y, z = cartesian_from_spherical([r, θ, φ])

    r_center = sqrt(l_center^2 + b0^2 - Q^2)

    out[1] = (x - r_center)^2 + y^2 + z^2 - star_radius^2
    out[2] = (rmax - r) * (rmax + r)
end

function star_across_wormhole_affect!(integrator, idx)
    if integrator.u[2] < 0.0 || idx == 2
        terminate!(integrator)
    end
end
