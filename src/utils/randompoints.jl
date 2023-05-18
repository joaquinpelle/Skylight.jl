function random_uniform_points_annulus!(points, rmin, rmax, coord_system)

    N = size(points,2)

    r = random_cylindrical_radius(N, rmin, rmax)
    φ = random_azimuthal_angle(N)

    set_points_on_equatorial_plane!(points, r, φ, coord_system)

end

function random_uniform_points_disk!(points, rmax, coord_system)

    random_uniform_points_annulus!(points, 0.0, rmax, coord_system)

end

function random_uniform_points_unit_spherical_cap!(points, θmax_in_degrees, coord_system)

    N = size(points,2)

    θmax_in_radians = deg2rad(θmax_in_degrees)

    θ = random_polar_angle(N, θmax_in_radians)
    φ = random_azimuthal_angle(N)
    
    set_points_on_unit_sphere!(points, θ, φ, coord_system)

end

function random_uniform_points_unit_sphere!(points, coord_system)

    random_uniform_points_unit_spherical_cap!(points, π, coord_system)

end

function random_uniform_points_unit_hemisphere!(points, coord_system)

    random_uniform_points_unit_spherical_cap!(points, π/2, coord_system)

end

function random_uniform_points_unit_hemisphere_xaxis!(v,::CartesianClass)

    random_uniform_points_unit_hemisphere!(v, CartesianClass())
    rotate_around_y_axis!(v,90)

end

random_cylindrical_radius(N, rmin, rmax) = sqrt.(rmin^2 .+ (rmax^2 - rmin^2)*rand(N))
random_polar_angle(N, θmax) = acos.(1.0.-(1.0-cos(θmax))*rand(N))
random_azimuthal_angle(N) = 2π*rand(N)

function set_points_on_equatorial_plane!(points, r, φ, ::CartesianClass)
    
    @. begin
    
        points[1,:] = r*cos(φ)
        points[2,:] = r*sin(φ)
    
    end

end

function set_points_on_equatorial_plane!(points, r, φ, ::SphericalClass)

    @. begin

        points[1,:] = r
        points[2,:] = π/2
        points[3,:] = φ        

    end

end

function set_points_on_unit_sphere!(points, θ, φ, ::CartesianClass)

    @. begin

        points[1,:] = sin(θ)*cos(φ)
        points[2,:] = sin(θ)*sin(φ)
        points[3,:] = cos(θ)  
    
    end

end

function set_points_on_unit_sphere!(points, θ, φ, ::SphericalClass)

    @. begin

        points[1,:] = 1.0
        points[2,:] = θ
        points[3,:] = φ        

    end

end


