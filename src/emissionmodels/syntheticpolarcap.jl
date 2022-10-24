export SyntheticPolarCap

@with_kw struct SyntheticPolarCap <: OpaqueInteriorSurfaceEmissionModel

    number_of_points::Int64
    star_radius::Float64
    angular_speed::Float64
    misalignment_angle_in_degrees::Float64 
    angular_radius_in_degrees::Float64
    temperature::Float64
    angular_radius_in_radians::Float64 = deg2rad(angular_radius_in_degrees)
    misalignment_angle_in_radians::Float64 = deg2rad(misalignment_angle_in_degrees)
    
end

function set_surface_differential!(covector, position, model::SyntheticPolarCap, coord_system::CartesianKind)

    @views begin
        t,x,y,z = position
        R = model.star_radius
    end

    covector[1] = 0.0
    covector[2] = 2x
    covector[3] = 2y
    covector[4] = 2z

end

function set_model_four_velocity!(vector, position, metric, model::SyntheticPolarCap, coord_system)
        
    angular_speed = model.angular_speed
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coord_system)
    
end

function get_space_positions(model::SyntheticPolarCap, coord_system::CartesianKind)

    space_positions = zeros(3, model.number_of_points)

    random_uniform_points_unit_spherical_cap!(space_positions, model.angular_radius_in_degrees, coord_system)

    rotate_around_y_axis!(space_positions, model.misalignment_angle_in_degrees)

    space_positions .*= model.star_radius

    return space_positions

end

get_number_of_points(model::SyntheticPolarCap) = model.number_of_points

get_callback(model::SyntheticPolarCap, coord_system::CartesianKind) = star_cartesian_callback()
get_callback(model::SyntheticPolarCap, coord_system::SphericalKind) = star_spherical_callback()