function random_uniform_points_annulus!(points, rmin, rmax, coords_top)
    N = size(points,2)
    r = random_cylindrical_radius(N, rmin, rmax)
    φ = random_azimuthal_angle(N)
    set_points_on_equatorial_plane!(points, r, φ, coords_top)
    return nothing
end

function random_uniform_points_disk!(points, rmax, coords_top)
    random_uniform_points_annulus!(points, 0.0, rmax, coords_top)
    return nothing
end

function random_uniform_points_unit_spherical_cap!(points, θmax_in_degrees, coords_top)
    N = size(points,2)
    θmax_in_radians = deg2rad(θmax_in_degrees)
    θ = random_polar_angle(N, θmax_in_radians)
    φ = random_azimuthal_angle(N)
    set_points_on_unit_sphere!(points, θ, φ, coords_top)
    return nothing
end

function random_uniform_points_unit_sphere!(points, coords_top)
    random_uniform_points_unit_spherical_cap!(points, π, coords_top)
    return nothing
end

function random_uniform_points_unit_hemisphere!(points, coords_top)
    random_uniform_points_unit_spherical_cap!(points, π/2, coords_top)
    return nothing
end

function random_uniform_points_unit_hemisphere_xaxis!(v,::CartesianTopology)
    random_uniform_points_unit_hemisphere!(v, CartesianTopology())
    rotate_around_y_axis!(v,90)
    return nothing
end

random_cylindrical_radius(N, rmin, rmax) = sqrt.(rmin^2 .+ (rmax^2 - rmin^2)*rand(N))
random_polar_angle(N, θmax) = acos.(1.0.-(1.0-cos(θmax))*rand(N))
random_azimuthal_angle(N) = 2π*rand(N)

function set_points_on_equatorial_plane!(points, r, φ, ::CartesianTopology)
    @. begin
        points[1,:] = r*cos(φ)
        points[2,:] = r*sin(φ)
    end
    return nothing
end

function set_points_on_equatorial_plane!(points, r, φ, ::SphericalTopology)
    @. begin
        points[1,:] = r
        points[2,:] = π/2
        points[3,:] = φ        
    end
    return nothing
end

function set_points_on_unit_sphere!(points, θ, φ, ::CartesianTopology)
    @. begin
        points[1,:] = sin(θ)*cos(φ)
        points[2,:] = sin(θ)*sin(φ)
        points[3,:] = cos(θ)  
    end
    return nothing
end

function set_points_on_unit_sphere!(points, θ, φ, ::SphericalTopology)
    @. begin
        points[1,:] = 1.0
        points[2,:] = θ
        points[3,:] = φ        
    end
    return nothing
end