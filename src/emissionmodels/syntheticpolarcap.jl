struct SyntheticPolarCap <: EmissionModel end

@with_kw struct SyntheticPolarCapParameters

    Npoints::Int64
    NS_radius::Float64
    misalignment_angle_in_degrees::Float64 
    angular_radius_in_degrees::Float64
    temperature::Float64
    angular_radius_in_radians::Float64 = deg2rad(angular_radius_in_degrees)
    misalignment_angle_in_radians::Float64 = deg2rad(misalignment_angle_in_degrees)
    
end

function synthetic_polar_cap(par::SyntheticPolarCapParameters)

    dataframe = zeros(4, par.Npoints)

    @views begin
        points = dataframe[1:3,:]
        temperatures = dataframe[4,:]
    end
    
    random_isotropic_unit_vectors_cone!(points, par.angular_radius_in_degrees)
    
    @. temperatures = par.temperature
    
    rotate_around_y_axis!(points, par.misalignment_angle_in_degrees)

    return dataframe

end
