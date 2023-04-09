export BosonStarAccretionDisk 

@with_kw struct BosonStarAccretionDisk{T} <: SurfaceEmissionModel

    inner_radius::Float64
    outer_radius::Float64

    temperature_file::String
    temperature_interpolator::T = build_temperature_interpolator(temperature_file)

end

function set_surface_differential!(covector, position, model::BosonStarAccretionDisk, coord_system::SphericalClass)

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0
  
end

function set_emitter_four_velocity!(vector, position, metric, spacetime, model::BosonStarAccretionDisk, coord_system)

    angular_speed = circular_geodesic_angular_speed(position, spacetime)
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coord_system)

end

function get_emitted_bolometric_flux(position, momentum, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, coord_system::SphericalClass)

    r = position[2]
    T = model.temperature_interpolator(r)
    
    return thermal_emission_bolometric_flux(T)

end

function get_emitted_specific_flux(position, momentum, energy, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, coord_system::SphericalClass)

    r = position[2]
    T = model.temperature_interpolator(r)
    
    return thermal_emission_specific_flux(energy, T)

end

function build_interpolator(temperature_file)

    data = readdlm(temperature_file, ' ', Float64, '\n')

    return my_interpolation(data[:,2], data[:,1], kind="cubicspline")

end


    
    

