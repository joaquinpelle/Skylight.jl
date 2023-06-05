contract(v, u) = dot(v, u)  
lower_index(v,metric) = metric*v
raise_index(v, metric_inverse) = metric_inverse*v
scalar_product(v,u,metric) = dot(v,metric,u)
norm_squared(v,metric) = scalar_product(v,v,metric)

"""Projection of v orthogonal to normalized timelike vector u"""
orthogonal_projection(v, u, metric) = v + u*scalar_product(v,u,metric) 

"""The cosine of the angle between vectors v and w as seen by observer at normalized four-velocity"""
function cos_angle_between_vectors(v, w, u, metric)
    vp = orthogonal_projection(v, u, metric)  
    wp = orthogonal_projection(w, u, metric)

    vp_wp = scalar_product(vp, wp, metric)
    vp2 = scalar_product(vp, vp, metric)
    wp2 = scalar_product(wp, wp, metric)

    return vp_wp/sqrt(vp2*wp2)
end

"""The cosine of the angle between null vectors v and w as seen by observer at normalized four-velocity u"""
cos_angle_between_null_vectors(v, w, u, metric) = 1.0+vector_scalar_product(v, w, metric)/(vector_scalar_product(v, u, metric)*vector_scalar_product(w, u, metric))

function normalize_timelike!(v, metric)
    v ./= sqrt(-norm_squared(v,metric))
    return nothing
end

function normalize_spacelike!(v, metric)
    v ./= sqrt(norm_squared(v,metric))
    return nothing
end

function spherical_from_cartesian(v)
    # Angles satisfy θ∈[0,π], φ∈[-π,π]
    r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    θ = acos(v[3]/r)      
    φ = atan(v[2],v[1])
    
    return [r,θ,φ]
end

function cartesian_from_spherical(v)
    x = v[1]*sin(v[2])*cos(v[3])
    y = v[1]*sin(v[2])*sin(v[3])
    z = v[1]*cos(v[2])
    
    return [x, y, z]
end

function rotate_around_y_axis!(v, angle_in_degrees)
    Nvectors = size(v, 2)

    ξ = deg2rad(angle_in_degrees)

    rotation_matrix = [cos(ξ) 0.0 sin(ξ); 0.0 1.0 0.0; -sin(ξ) 0.0 cos(ξ)]

    for i in 1:Nvectors

        v[:,i] .= rotation_matrix*v[:,i] 
    
    end
    return nothing
end