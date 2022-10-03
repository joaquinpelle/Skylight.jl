function random_uniform_points_disk!(points,radius)

    N = size(points,2)

    ρ = radius*sqrt.(rand(N))
    
    φ = 2π*rand(N)

    points[1,:] .= ρ.*cos.(φ)
    points[2,:] .= ρ.*sin.(φ)

    return points

end

function random_uniform_points_annulus!(points,inner_radius,outer_radius)

    N = size(points,2)

    ρ = sqrt.(inner_radius^2 .+ (outer_radius^2 - inner_radius^2)*rand(N))
    
    φ = 2π*rand(N)

    points[1,:] .= ρ.*cos.(φ)
    points[2,:] .= ρ.*sin.(φ)

    return points

end

function random_isotropic_unit_vectors_sphere!(v)

    N = size(v,2)
    
    θ = acos.(1.0.-2*rand(N))
    φ = 2π*(0.5.-rand(N))
    
    @. begin

        v[1,:] = sin(θ)*cos(φ)
        v[2,:] = sin(θ)*sin(φ)
        v[3,:] = cos(θ)
    
    end

    return v

end

function random_isotropic_unit_vectors_hemisphere!(v)

    N = size(v,2)
    
    θ = acos.(1.0.-rand(N))
    φ = 2π*(0.5.-rand(N))
    
    @. begin

        v[1,:] = sin(θ)*cos(φ)
        v[2,:] = sin(θ)*sin(φ)
        v[3,:] = cos(θ)
    
    end

    return v

end


function random_isotropic_unit_vectors_cone!(v,angular_radius_in_degrees)

    N = size(v,2)

    α = deg2rad(angular_radius_in_degrees)

    θ = acos.(1.0.-(1.0-cos(α))*rand(N))
    φ = 2π*(0.5.-rand(N))
    
    @. begin

        v[1,:] = sin(θ)*cos(φ)
        v[2,:] = sin(θ)*sin(φ)
        v[3,:] = cos(θ)
    
    end

    return v

end
