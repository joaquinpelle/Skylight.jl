contract(v, u) = dot(v, u)
lower_index(v, metric) = metric * v
raise_index(v, metric_inverse) = metric_inverse * v
scalar_product(v, u, metric) = dot(v, metric, u)
norm_squared(v, metric) = scalar_product(v, v, metric)

function normalize!(v, metric)
    v ./= sqrt(abs(norm_squared(v, metric)))
    return nothing
end

function normalize_timelike!(v, metric)
    v ./= sqrt(-norm_squared(v, metric))
    return nothing
end

function normalize_spacelike!(v, metric)
    v ./= sqrt(norm_squared(v, metric))
    return nothing
end

is_timelike(v, metric) = norm_squared(v, metric) < 0.0
is_spacelike(v, metric) = norm_squared(v, metric) > 0.0

"""Projection of v orthogonal to normalized timelike vector u"""
orthogonal_projection(v, u, metric) = v + u * scalar_product(v, u, metric)

"""The cosine of the angle between vectors v and w as seen by observer at normalized four-velocity"""
function cos_angle_between_vectors(v, w, u, metric)
    vp = orthogonal_projection(v, u, metric)
    wp = orthogonal_projection(w, u, metric)
    vp_wp = scalar_product(vp, wp, metric)
    vp2 = scalar_product(vp, vp, metric)
    wp2 = scalar_product(wp, wp, metric)
    return vp_wp / sqrt(vp2 * wp2)
end

"""The cosine of the angle between null vectors v and w as seen by observer at normalized four-velocity u"""
function cos_angle_between_null_vectors(v, w, u, metric)
    1.0 +
    vector_scalar_product(v, w, metric) /
    (vector_scalar_product(v, u, metric) * vector_scalar_product(w, u, metric))
end

function spherical_from_cartesian(v)
    # Angles satisfy θ∈[0,π], φ∈[-π,π]
    r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    θ = acos(v[3] / r)
    φ = atan(v[2], v[1])
    return SVector{3}(r, θ, φ)
end

function cartesian_from_spherical(v)
    x = v[1] * sin(v[2]) * cos(v[3])
    y = v[1] * sin(v[2]) * sin(v[3])
    z = v[1] * cos(v[2])
    return SVector{3}(x, y, z)
end

function rotate_around_y_axis!(v, angle_in_degrees)
    ξ = deg2rad(angle_in_degrees)
    rotation_matrix = [cos(ξ) 0.0 sin(ξ); 0.0 1.0 0.0; -sin(ξ) 0.0 cos(ξ)]
    for i in axes(v, 2)
        v[:, i] .= rotation_matrix * v[:, i]
    end
    return nothing
end

function has_lorentzian_signature(metric)
    is_timelike(SVector{4}(1.0, 0.0, 0.0, 0.0), metric) &&
        is_spacelike(SVector{4}(0.0, 1.0, 0.0, 0.0), metric) &&
        is_spacelike(SVector{4}(0.0, 0.0, 1.0, 0.0), metric) &&
        is_spacelike(SVector{4}(0.0, 0.0, 0.0, 1.0), metric)
end

function equatorial_position(r, ::SphericalTopology)
    return SVector{4}(0.0, r, π / 2, 0.0)
end

function equatorial_position(r::Real, φ::Real, ::SphericalTopology)
    return SVector{4}(0.0, r, π / 2, φ)
end

function equatorial_position(r, ::CartesianTopology)
    spherical_position = equatorial_position(r, SphericalTopology())
    return SVector{4}(0.0, cartesian_from_spherical(spherical_position[2:end])...)
end

function equatorial_position(r::Real, φ::Real, ::CartesianTopology)
    spherical_position = equatorial_position(r, φ, SphericalTopology())
    return SVector{4}(0.0, cartesian_from_spherical(spherical_position[2:end])...)
end

function equatorial_position!(position, r, ::SphericalTopology)
    position[2] = r
    position[3] = π / 2
    return nothing
end

function equatorial_position!(position, r, ::CartesianTopology)
    position[2] = r
    position[3] = 0.0
    position[4] = 0.0
    return nothing
end
