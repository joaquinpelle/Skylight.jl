contract(vector, covector)  = vector'covector  
lower_index(vector, metric) = metric*vector 
scalar_product(v,u,metric) = v'*metric*u
norm_squared(v,metric) = scalar_product(v,v,metric)

∂t() = [1.0, 0.0, 0.0, 0.0]

function cartesian_to_spherical(v)

    # Angles satisfy θ∈[0,π], φ∈[-π,π]

    r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    θ = acos(v[3]/r)      
    φ = atan(v[2],v[1])
    
    return [r,θ,φ]
    
end