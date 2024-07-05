"""
    CircularHotSpot <: AbstractSurfaceEmissionModel

Circular hot spot on the surface of a neutron star.

# Fields
- `star_radius_in_km::Float64`: The radius of the star in km. Must be positive.
- `spin_frequency_in_Hz::Float64`: The spin frequency of the star in Hz. Must be non-zero.
- `center_colatitude_in_degrees::Float64`: The colatitude of the spot center in degrees. Must be in the range [0, 90].
- `angular_radius_in_radians::Float64`: The angular radius of the polar cap, in radians. Must be positive.
- `temperature_in_keV::Float64`: The temperature of the polar cap. Must be positive.

# Examples
```julia
hot_spot = CircularHotSpot(
    star_radius_in_km = 12.0,
    spin_frequency_in_Hz = 200,
    center_colatitude_in_degrees = 30.0,
    angular_radius_in_radians = 1.0,
    temperature_in_keV = 0.35
    M1 = 1.4
)
```
"""
@with_kw struct CircularHotSpot <: AbstractSurfaceEmissionModel
    star_radius_in_km::Float64
    spin_frequency_in_Hz::Float64
    center_colatitude_in_degrees::Float64
    angular_radius_in_radians::Float64
    temperature_in_keV::Float64
    M1::Float64 
    star_radius::Float64 = CGS_to_geometrized(1e5*star_radius_in_km, Dimensions.length; M1 = M1)
    angular_speed::Float64 = CGS_to_geometrized(2π*spin_frequency_in_Hz, Dimensions.frequency; M1 = M1)
    angular_radius_in_degrees::Float64 = rad2deg(angular_radius_in_radians)
    center_colatitude_in_radians::Float64 = deg2rad(center_colatitude_in_degrees)
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

function system_period(model::CircularHotSpot)
    return 2π / model.angular_speed
end

function surface_differential!(covector, position, ::CircularHotSpot, ::SphericalTopology)
    covector[1] = 0.0
    covector[2] = 1.0
    covector[3] = 0.0
    covector[4] = 0.0
    return nothing
end

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

function rest_frame_specific_intensity(position, 
    momentum, 
    energy, 
    rest_frame_four_velocity, 
    metric, 
    spacetime, 
    model::CircularHotSpot, 
    coords_top)
    return thermal_emission_specific_intensity(energy, model.temperature)
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

function space_positions(npoints, spacetime, model::CircularHotSpot, ::SphericalTopology, cache)
    space_pos = space_positions(npoints, spacetime, model, CartesianTopology(), cache)
    for v in eachcol(space_pos)
        v .= spherical_from_cartesian(v)
    end
    return space_pos
end

function photon_package_weight(position, 
    momentum, 
    emitted_energy, 
    metric,
    rest_frame_four_velocity,
    surface_normal, 
    spacetime,
    model::CircularHotSpot, 
    coords_top)
    Iem = rest_frame_specific_intensity(position, momentum, emitted_energy, rest_frame_four_velocity, metric, spacetime, model, coords_top)
    cosθem = cos_angle_between_vectors(momentum, surface_normal, rest_frame_four_velocity, metric)
    return Iem*cosθem/emitted_energy
end