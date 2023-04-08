export BosonStarAccretionDisk 

@with_kw struct BosonStarAccretionDisk <: SurfaceEmissionModel

    inner_radius::Float64
    outer_radius::Float64

    temperature_file::String

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

function get_emitted_bolometric_flux(position, momentum, emitter_four_velocity, metric, spacetime, model::BosonStarAccretionDisk, coord_system)

    T = get_temperature(position, model, coord_system)

    return thermal_emission_bolometric_flux(T)

end


function get_temperature(position, model::BosonStarAccretionDisk, coord_system::SphericalClass)

    t, r, θ, φ = position

    return model.temperature_interpolator(r)

end

function build_interpolator(temperature)



end

function temperature(position, model::BosonStarAccretionDisk)






function read_temperature_file(model::BosonStarAccretionDisk)

    temperature = readdlm(model.temperature_file, ' ', Float64, '\n')

end



    
    

