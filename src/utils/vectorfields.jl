∂t() = [1.0, 0.0, 0.0, 0.0]
∂φ() = [0.0, 0.0, 0.0, 1.0]

time_translation_generator() = ∂t()

function rotation_generators(position, ::CartesianTopology)
    t, x, y, z = position
    X = [0.0, 0.0, z, -y]
    Y = [0.0,  z, 0.0,-x]
    Z = [0.0,  y, -x, 0.0]
    return [X, Y, Z]
end

function zaxis_rotation_generator(position, ::CartesianTopology)
    t, x, y, z = position
    return [0.0, y, -x, 0.0]
end

function rotation_generators(position, ::SphericalTopology)
    t, r, θ, φ = position
    X = [0.0, 0.0, sin(φ), cot(θ)*cos(φ)]
    Y = [0.0, 0.0, cos(φ),-cot(θ)*sin(φ)]
    return [X, Y, ∂φ()]
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
    
    vectors[3,1] = z
    vectors[4,1] = -y

    vectors[2,2] = z
    vectors[4,2] = -x

    vectors[2,3] = y
    vectors[3,3] = -x

    return nothing
end

function zaxis_rotation_generator!(v, position, ::CartesianTopology)
    t, x, y, z = position
    v[2] = y
    v[3] = -x
    return nothing
end

function rotation_generators!(vectors, position, ::SphericalTopology)
    fill!(vectors, 0.0)
    t, r, θ, φ = position
    vectors[3,1] = sin(φ)
    vectors[4,1] = cot(θ)*cos(φ)

    vectors[3,2] = cos(φ)
    vectors[4,2] = -cot(θ)*sin(φ)

    vectors[4,3] = 1.0
    return nothing 
end

function zaxis_rotation_generator!(v, position, ::SphericalTopology)
    fill!(v, 0.0)
    v[4] = 1.0
    return nothing
end


function tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, ::CartesianTopology)
    vector[1] =  1.0
    vector[2] = -angular_speed*position[3]
    vector[3] =  angular_speed*position[2]
    vector[4] =  0.0

    normalize_timelike!(vector,metric)
end

function tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, ::SphericalTopology)
    vector[1] =  1.0
    vector[2] =  0.0
    vector[3] =  0.0
    vector[4] =  angular_speed

    normalize_timelike!(vector,metric)
end
