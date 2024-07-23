function callback_parameters(spacetime::AbstractBlackHoleSpacetime,
    ::DexterAgolDisk,
    configurations;
    rhorizon_bound)
    rhorizon = event_horizon_radius(spacetime)
    rmin = rhorizon + rhorizon_bound
    rmax = max_radius(configurations)
    return BlackHoleEquatorialPlaneCallbackParameters(rmin = rmin,
        rmax = rmax)
end

@with_kw struct BlackHoleEquatorialPlaneCallbackParameters <: AbstractCallbackParameters
    rmin::Float64
    rmax::Float64
    @assert rmin<rmax "rmin must be smaller than rmax"
end

function callback(::AbstractBlackHoleSpacetime,
    ::DexterAgolDisk,
    ::SphericalTopology)
    black_hole_equatorial_plane_spherical_coordinates_callback()
end
function callback(::SchwarzschildSpacetimeKerrSchildCoordinates,
    ::DexterAgolDisk,
    ::CartesianTopology)
    nonspinning_black_hole_equatorial_plane_cartesian_coordinates_callback()
end
function callback(::KerrSpacetimeKerrSchildCoordinates,
    ::DexterAgolDisk,
    ::CartesianTopology)
    spinning_black_hole_equatorial_plane_cartesian_coordinates_callback()
end

function spinning_black_hole_equatorial_plane_cartesian_coordinates_callback()
    VectorContinuousCallback(spinning_black_hole_equatorial_plane_cartesian_coordinates_condition,
        black_hole_equatorial_plane_affect!,
        2)
end
function nonspinning_black_hole_equatorial_plane_cartesian_coordinates_callback()
    VectorContinuousCallback(nonspinning_black_hole_equatorial_plane_cartesian_coordinates_condition,
        black_hole_equatorial_plane_affect!,
        2)
end
function black_hole_equatorial_plane_spherical_coordinates_callback()
    VectorContinuousCallback(black_hole_equatorial_plane_spherical_coordinates_condition,
        black_hole_equatorial_plane_affect!,
        2)
end

function black_hole_equatorial_plane_affect!(integrator,
    idx)
    terminate!(integrator)
end

function nonspinning_black_hole_equatorial_plane_cartesian_coordinates_condition(out,
    u,
    t,
    integrator)
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    r2 = u[2] * u[2] + u[3] * u[3] + u[4] * u[4]
    out[1] = u[4]
    out[2] = (rmax * rmax - r2) * (r2 - rmin * rmin)
end

function spinning_black_hole_equatorial_plane_cartesian_coordinates_condition(out,
    u,
    t,
    integrator)
    a = integrator.p.spacetime.a
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    rho2_a2 = u[2] * u[2] + u[3] * u[3] + u[4] * u[4] - a * a
    r2 = 0.5 * (rho2_a2 + sqrt(rho2_a2 * rho2_a2 + 4.0 * a * a * u[4] * u[4]))
    out[1] = u[4]
    out[2] = (rmax * rmax - r2) * (r2 - rmin * rmin)
end

function black_hole_equatorial_plane_spherical_coordinates_condition(out, u, t, integrator)
    rmax = integrator.p.cbp.rmax
    rmin = integrator.p.cbp.rmin
    out[1] = u[3] - π / 2
    out[2] = (rmax - u[2]) * (u[2] - rmin)
end