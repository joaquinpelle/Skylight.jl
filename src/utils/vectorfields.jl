∂t() = [1.0, 0.0, 0.0, 0.0]
∂φ() = [0.0, 0.0, 0.0, 1.0]

#Generators
time_translation_generator() = ∂t()

function rotation_generators(position, ::CartesianTopology)
    t, x, y, z = position
    X = SVector{4}(0.0, 0.0, z, -y)
    Y = SVector{4}(0.0, z, 0.0, -x)
    Z = SVector{4}(0.0, y, -x, 0.0)
    return (X, Y, Z)
end

function zaxis_rotation_generator(position, ::CartesianTopology)
    t, x, y, z = position
    return SVector{4}(0.0, y, -x, 0.0)
end

function rotation_generators(position, ::SphericalTopology)
    t, r, θ, φ = position
    X = SVector{4}(0.0, 0.0, sin(φ), cot(θ) * cos(φ))
    Y = SVector{4}(0.0, 0.0, cos(φ), -cot(θ) * sin(φ))
    return (X, Y, ∂φ())
end

zaxis_rotation_generator(position, ::SphericalTopology) = ∂φ()

function time_translation_generator!(v)
    v[1] = 1.0
    v[2] = 0.0
    v[3] = 0.0
    v[4] = 0.0
    return nothing
end

function rotation_generators!(vectors, position, ::CartesianTopology)
    fill!(vectors, 0.0)
    t, x, y, z = position

    vectors[3, 1] = z
    vectors[4, 1] = -y

    vectors[2, 2] = z
    vectors[4, 2] = -x

    vectors[2, 3] = y
    vectors[3, 3] = -x
    return nothing
end

function zaxis_rotation_generator!(v, position, ::CartesianTopology)
    fill!(v, 0.0)
    t, x, y, z = position
    v[2] = y
    v[3] = -x
    return nothing
end

function rotation_generators!(vectors, position, ::SphericalTopology)
    fill!(vectors, 0.0)
    t, r, θ, φ = position
    vectors[3, 1] = sin(φ)
    vectors[4, 1] = cot(θ) * cos(φ)

    vectors[3, 2] = cos(φ)
    vectors[4, 2] = -cot(θ) * sin(φ)

    vectors[4, 3] = 1.0
    return nothing
end

function zaxis_rotation_generator!(v, position, ::SphericalTopology)
    fill!(v, 0.0)
    v[4] = 1.0
    return nothing
end

#Four-velocities
function static_four_velocity(metric)
    vector = time_translation_generator()
    normalize_timelike!(vector, metric)
    return vector
end

function static_four_velocity!(vector, metric)
    time_translation_generator!(vector)
    normalize_timelike!(vector, metric)
    return nothing
end

function static_four_velocity_allowing_spacelike!(vector, metric)
    time_translation_generator!(vector)
    normalize!(vector, metric)
    return nothing
end

function circular_motion_four_velocity(position, angular_speed, metric, ::CartesianTopology)
    vector = zeros(4)
    vector[1] = 1.0
    vector[2] = -angular_speed * position[3]
    vector[3] = angular_speed * position[2]
    normalize_timelike!(vector, metric)
    return nothing
end

function circular_motion_four_velocity(position, angular_speed, metric, ::SphericalTopology)
    vector = zeros(4)
    vector[1] = 1.0
    vector[4] = angular_speed
    normalize_timelike!(vector, metric)
    return nothing
end

function circular_motion_four_velocity!(vector,
    position,
    angular_speed,
    metric,
    ::CartesianTopology)
    vector[1] = 1.0
    vector[2] = -angular_speed * position[3]
    vector[3] = angular_speed * position[2]
    vector[4] = 0.0
    normalize_timelike!(vector, metric)
    return nothing
end

function circular_motion_four_velocity!(vector,
    position,
    angular_speed,
    metric,
    ::SphericalTopology)
    vector[1] = 1.0
    vector[2] = 0.0
    vector[3] = 0.0
    vector[4] = angular_speed
    normalize_timelike!(vector, metric)
    return nothing
end

function circular_motion_four_velocity_allowing_spacelike!(vector,
    position,
    angular_speed,
    metric,
    ::CartesianTopology)
    vector[1] = 1.0
    vector[2] = -angular_speed * position[3]
    vector[3] = angular_speed * position[2]
    vector[4] = 0.0
    normalize!(vector, metric)
    return nothing
end

function circular_motion_four_velocity_allowing_spacelike!(vector,
    position,
    angular_speed,
    metric,
    ::SphericalTopology)
    vector[1] = 1.0
    vector[2] = 0.0
    vector[3] = 0.0
    vector[4] = angular_speed
    normalize!(vector, metric)
    return nothing
end
