@with_kw struct SyntheticPolarCap <: AbstractSurfaceEmissionModel
    star_radius::Float64
    angular_speed::Float64
    misalignment_angle_in_degrees::Float64
    angular_radius_in_degrees::Float64
    temperature::Float64
    angular_radius_in_radians::Float64 = deg2rad(angular_radius_in_degrees)
    misalignment_angle_in_radians::Float64 = deg2rad(misalignment_angle_in_degrees)

    @assert star_radius>0.0 "star_radius must be positive"
    @assert angular_speed!=0.0 "angular_speed must be non-zero"
    @assert angular_radius_in_degrees>0.0 "angular_radius_in_degrees must be positive"
    @assert misalignment_angle_in_degrees>=0.0 "misalignment_angle_in_degrees must be non-negative"
    @assert misalignment_angle_in_degrees<=90.0 "misalignment_angle_in_degrees must be smaller than 90.0"
    @assert temperature>0.0 "temperature must be positive"
end

opaque_interior_surface_trait(::SyntheticPolarCap) = IsOpaqueInteriorSurface()
stationarity(::SyntheticPolarCap) = IsStationary()

function surface_differential!(covector, position, ::SyntheticPolarCap, ::CartesianTopology)
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
    model::SyntheticPolarCap,
    coords_top)
    angular_speed = model.angular_speed
    circular_motion_four_velocity!(vector, position, angular_speed, metric, coords_top)
    return nothing
end

function space_positions(npoints, spacetime, model::SyntheticPolarCap, coords_top::CartesianTopology, cache)
    space_pos = zeros(3, npoints)
    random_uniform_points_unit_spherical_cap!(space_pos,
        model.angular_radius_in_degrees,
        coords_top)
    rotate_around_y_axis!(space_pos, model.misalignment_angle_in_degrees)
    space_pos .*= model.star_radius
    return space_pos
end
