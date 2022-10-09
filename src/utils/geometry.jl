contract(v, u) = dot(v, u)  
lower_index(v,metric) = metric*v
raise_index(v, metric_inverse) = metric_inverse*v
scalar_product(v,u,metric) = dot(v,metric,u)
norm_squared(v,metric) = scalar_product(v,v,metric)

∂t() = [1.0, 0.0, 0.0, 0.0]

function normalize_timelike!(v, metric)

    v ./= sqrt(-norm_squared(v,metric))

end

function normalize_spacelike!(v, metric)

    v ./= sqrt(norm_squared(v,metric))

end

function tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coord_system::CartesianKind)
    
    vector[1] =  1.0
    vector[2] = -angular_speed*position[3]
    vector[3] =  angular_speed*position[2]
    vector[4] =  0.0

    normalize_timelike!(vector,metric)

end

function tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coord_system::SphericalKind)
    
    vector[1] =  1.0
    vector[2] =  0.0
    vector[3] =  0.0
    vector[4] =  angular_speed

    normalize_timelike!(vector,metric)

end

function spherical_from_cartesian(v)

    # Angles satisfy θ∈[0,π], φ∈[-π,π]

    r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    θ = acos(v[3]/r)      
    φ = atan(v[2],v[1])
    
    return [r,θ,φ]
    
end

function rotate_around_y_axis!(v, angle_in_degrees)

    Nvectors = size(v, 2)

    ξ = deg2rad(angle_in_degrees)

    rotation_matrix = [cos(ξ) 0.0 sin(ξ); 0.0 1.0 0.0; -sin(ξ) 0.0 cos(ξ)]

    for i in 1:Nvectors

        v[:,i] .= rotation_matrix*v[:,i] 
    
    end

end