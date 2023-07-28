"""
    random_uniform_points_annulus!(points, rmin, rmax, coords_top)

Populate the given array `points` with points randomly and uniformly 
distributed in an annulus with inner radius `rmin` and outer radius `rmax`.
The type of the coordinate system is determined by `coords_top`.
"""
function random_uniform_points_annulus!(points, rmin, rmax, coords_top)
    N = size(points,2)
    r = random_cylindrical_radius(N, rmin, rmax)
    φ = random_azimuthal_angle(N)
    points_on_equatorial_plane!(points, r, φ, coords_top)
    return nothing
end

"""
    random_uniform_points_disk!(points, rmax, coords_top)

Populate the given array `points` with points randomly and uniformly 
distributed in a disk with radius `rmax`.
The type of the coordinate system is determined by `coords_top`.
"""
function random_uniform_points_disk!(points, rmax, coords_top)
    random_uniform_points_annulus!(points, 0.0, rmax, coords_top)
    return nothing
end

"""
    random_uniform_points_unit_spherical_cap!(points, θmax_in_degrees, coords_top)

Populate the given array `points` with points randomly and uniformly 
distributed in a spherical cap with maximum polar angle `θmax_in_degrees` degrees.
The type of the coordinate system is determined by `coords_top`.
"""
function random_uniform_points_unit_spherical_cap!(points, θmax_in_degrees, coords_top)
    N = size(points,2)
    θmax_in_radians = deg2rad(θmax_in_degrees)
    θ = random_polar_angle(N, θmax_in_radians)
    φ = random_azimuthal_angle(N)
    points_on_unit_sphere!(points, θ, φ, coords_top)
    return nothing
end

"""
    random_uniform_points_unit_sphere!(points, coords_top)

Populate the given array `points` with points randomly and uniformly 
distributed on the unit sphere.
The type of the coordinate system is determined by `coords_top`.
"""
function random_uniform_points_unit_sphere!(points, coords_top)
    random_uniform_points_unit_spherical_cap!(points, 180, coords_top)
    return nothing
end

"""
    random_uniform_points_unit_hemisphere!(points, coords_top)

Populate the given array `points` with points randomly and uniformly 
distributed on the unit hemisphere.
The type of the coordinate system is determined by `coords_top`.
"""
function random_uniform_points_unit_hemisphere!(points, coords_top)
    random_uniform_points_unit_spherical_cap!(points, 90, coords_top)
    return nothing
end

"""
    random_uniform_points_unit_hemisphere_xaxis!(v,::CartesianTopology)

Populate the given array `v` with points randomly and uniformly 
distributed on the unit hemisphere along the x-axis.
"""
function random_uniform_points_unit_hemisphere_xaxis!(v,::CartesianTopology)
    random_uniform_points_unit_hemisphere!(v, CartesianTopology())
    rotate_around_y_axis!(v,90)
    return nothing
end

"""
    random_cylindrical_radius(N, rmin, rmax)

Generates `N` random numbers representing cylindrical radii, following a 
uniform distribution between `rmin` and `rmax`.
"""
random_cylindrical_radius(N, rmin, rmax) = sqrt.(rmin^2 .+ (rmax^2 - rmin^2)*rand(N))

"""
    random_polar_angle(N, θmax)

Generates `N` random numbers representing polar angles, following a 
uniform distribution from 0 to `θmax` (in radians).
"""
random_polar_angle(N, θmax) = acos.(1.0.-(1.0-cos(θmax))*rand(N))

"""
    random_azimuthal_angle(N)

Generates `N` random numbers representing azimuthal angles, following a 
uniform distribution from 0 to 2π.
"""
random_azimuthal_angle(N) = 2π*rand(N)

"""
    points_on_equatorial_plane!(points, r, φ, ::CartesianTopology)

Sets the points on the equatorial plane from spherical coordinates
"""
function points_on_equatorial_plane!(points, r, φ, ::CartesianTopology)
    @. begin
        points[1,:] = r*cos(φ)
        points[2,:] = r*sin(φ)
    end
    return nothing
end

function points_on_equatorial_plane!(points, r, φ, ::SphericalTopology)
    @. begin
        points[1,:] = r
        points[2,:] = π/2
        points[3,:] = φ        
    end
    return nothing
end

"""
    points_on_unit_sphere!(points, θ, φ, ::AbstractCoordinatesTopology)

Transforms the points on the unit sphere from the angular coordinates
"""
function points_on_unit_sphere!(points, θ, φ, ::CartesianTopology)
    @. begin
        points[1,:] = sin(θ)*cos(φ)
        points[2,:] = sin(θ)*sin(φ)
        points[3,:] = cos(θ)  
    end
    return nothing
end

function points_on_unit_sphere!(points, θ, φ, ::SphericalTopology)
    @. begin
        points[1,:] = 1.0
        points[2,:] = θ
        points[3,:] = φ        
    end
    return nothing
end