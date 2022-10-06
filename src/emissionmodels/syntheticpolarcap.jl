@with_kw struct SyntheticPolarCap <: EmissionModel

    number_of_points::Int64
    NS_radius::Float64
    angular_speed::Float64
    misalignment_angle_in_degrees::Float64 
    angular_radius_in_degrees::Float64
    temperature::Float64
    angular_radius_in_radians::Float64 = deg2rad(angular_radius_in_degrees)
    misalignment_angle_in_radians::Float64 = deg2rad(misalignment_angle_in_degrees)
    
end

get_number_of_points(model::SyntheticPolarCap) = model.number_of_points

function set_local_four_velocity!(vector, position, gμν, model::SyntheticPolarCap)

    angular_speed = model.angular_speed
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, gμν, CartesianKind())
    
end

function synthetic_polar_cap_dataframe(model::SyntheticPolarCap)

    dataframe = zeros(4, model.number_of_points)

    @views begin
        points = dataframe[1:3,:]
        temperatures = dataframe[4,:]
    end
    
    random_uniform_points_unit_spherical_cap!(points, model.angular_radius_in_degrees, CartesianKind())
    
    @. temperatures = model.temperature
    
    rotate_around_y_axis!(points, model.misalignment_angle_in_degrees)

    return dataframe

end


