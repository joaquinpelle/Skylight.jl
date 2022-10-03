contract(vector, covector)  = vector'covector  
lower_index(vector, metric) = metric*vector 
scalar_product(v,u,metric) = v'*metric*u
norm_squared(v,metric) = scalar_product(v,v,metric)

∂t() = [1.0, 0.0, 0.0, 0.0]

function normalize_timelike!(v, metric)

    v ./= sqrt(-norm_squared(v,metric))

end

function normalize_spacelike!(v, metric)

    v ./= sqrt(norm_squared(v,metric))

end

function spherical_from_cartesian(v)

    # Angles satisfy θ∈[0,π], φ∈[-π,π]

    r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    θ = acos(v[3]/r)      
    φ = atan(v[2],v[1])
    
    return [r,θ,φ]
    
end

function rotate_around_y_axis!(v, angle_in_degrees)

    ξ = deg2rad(angle_in_degrees)

    rotation_matrix = [cos(ξ) 0.0 sin(ξ); 0.0 1.0 0.0; -sin(ξ) 0.0 cos(ξ)]

    for i in range(1, size(v,2))

        v[:,i] .= rotation_matrix*v[:,i] 
    
    end

end