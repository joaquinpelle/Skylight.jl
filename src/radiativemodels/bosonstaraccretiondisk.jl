@with_kw struct BosonStarAccretionDisk{T} <: AbstractSurfaceEmissionModel

    inner_radius::Float64
    outer_radius::Float64

    temperature_file::String
    temperature_interpolator::T = build_interpolator(temperature_file)

end

function set_emitter_four_velocity!(vector, position, metric, spacetime::BosonStarSpacetime, ::BosonStarAccretionDisk, coords_top)

    angular_speed = circular_geodesic_angular_speed(position, spacetime)
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coords_top)

end

function get_emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, ::SphericalTopology)

    r = position[2]
    T = model.temperature_interpolator(r)
    
    return thermal_emission_bolometric_intensity(T)

end

function get_emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, ::SphericalTopology)

    r = position[2]
    T = model.temperature_interpolator(r)
    
    return thermal_emission_specific_intensity(energy, T)

end

## The following function is used to check if the ray is inside the accretion disk
function is_final_position_at_source(position, spacetime, model::BosonStarAccretionDisk)

    r = position[2]

    return (r >= model.inner_radius) && (r <= model.outer_radius)
    
end    
    

