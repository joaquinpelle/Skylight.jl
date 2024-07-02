"""
    CircularHotSpot <: AbstractSurfaceEmissionModel

Circular hot spot on the surface of a neutron star.

# Fields
- `star_radius_in_km::Float64`: The radius of the star in km. Must be positive.
- `angular_speed_in_Hz::Float64`: The angular speed of the star in Hz. Must be non-zero.
- `center_colatitude_in_degrees::Float64`: The colatitude of the spot center in degrees. Must be in the range [0, 90].
- `angular_radius_in_radians::Float64`: The angular radius of the polar cap, in radians. Must be positive.
- `temperature::Float64`: The temperature of the polar cap. Must be positive.

# Examples
```julia
hot_spot = CircularHotSpot(
    star_radius_in_km = 12.0,
    angular_speed_in_Hz = 200,
    center_colatitude_in_degrees = 30.0,
    angular_radius_in_radians = 1.0,
    temperature_in_keV = 0.35
    M1 = 1.4
)
```
"""
@with_kw struct CircularHotSpot <: AbstractSurfaceEmissionModel
    star_radius_in_km::Float64
    angular_speed_in_Hz::Float64
    center_colatitude_in_degrees::Float64
    angular_radius_in_radians::Float64
    temperature_in_keV::Float64
    M1::Float64 
    star_radius::Float64 = CGS_to_geometrized(1e5*star_radius_in_km, Dimensions.length; M1 = M1)
    angular_speed::Float64 = CGS_to_geometrized(angular_speed_in_Hz, Dimensions.frequency; M1 = M1)
    angular_radius_in_degrees::Float64 = rad2deg(angular_radius_in_radians)
    center_colatitude_in_radians::Float64 = deg2rad(center_colatitude_in_radians)
    temperature::Float64 = keV_to_K(temperature_in_keV)

    @assert star_radius>0.0 "star_radius must be positive"
    @assert angular_speed!=0.0 "angular_speed must be non-zero"
    @assert angular_radius_in_radians>0.0 "angular_radius_in_radians must be positive"
    @assert center_colatitude_in_degrees>=0.0 "center_colatitude_in_degrees must be non-negative"
    @assert center_colatitude_in_degrees<=90.0 "center_colatitude_in_degrees must be smaller than 90.0"
    @assert temperature>0.0 "temperature must be positive"
end

opaque_interior_surface_trait(::CircularHotSpot) = IsOpaqueInteriorSurface()
stationarity(::CircularHotSpot) = IsStationary()

function surface_differential!(covector, position, ::CircularHotSpot, ::CartesianTopology)
    @views begin
        x = position[2]
        y = position[3]
        z = position[4]    
    end

    covector[1] = 0.0
    covector[2] = 2x
    covector[3] = 2y
    covector[4] = 2z
    return nothing
end

function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    model::CircularHotSpot,
    coords_top)
    angular_speed = model.angular_speed
    circular_motion_four_velocity!(vector, position, angular_speed, metric, coords_top)
    return nothing
end

function space_positions(npoints, spacetime, model::CircularHotSpot, coords_top::CartesianTopology, cache)
    space_pos = zeros(3, npoints)
    random_uniform_points_unit_spherical_cap!(space_pos,
        model.angular_radius_in_degrees,
        coords_top)
    rotate_around_y_axis!(space_pos, model.center_colatitude_in_degrees)
    space_pos .*= model.star_radius
    return space_pos
end