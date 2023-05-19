@with_kw struct BosonStarAccretionDisk{T} <: SurfaceEmissionModel

    inner_radius::Float64
    outer_radius::Float64

    temperature_file::String
    temperature_interpolator::NoSaveField{T} = NoSaveField(build_interpolator(temperature_file))

end

getproperty(model::BosonStarAccretionDisk, field) = getproperty_nosave(model, field)

function set_emitter_four_velocity!(vector, position, metric, spacetime::BosonStarSpacetime, ::BosonStarAccretionDisk, coord_system)

    angular_speed = circular_geodesic_angular_speed(position, spacetime)
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coord_system)

end

function get_emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, ::SphericalClass)

    r = position[2]
    T = model.temperature_interpolator(r)
    
    return thermal_emission_bolometric_intensity(T)

end

function get_emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, ::SphericalClass)

    r = position[2]
    T = model.temperature_interpolator(r)
    
    return thermal_emission_specific_intensity(energy, T)

end

## The following function is used to check if the ray is inside the accretion disk
function is_final_position_at_source(position, spacetime, model::BosonStarAccretionDisk)

    r = position[2]

    return (r >= model.inner_radius) && (r <= model.outer_radius)
    
end    
    

