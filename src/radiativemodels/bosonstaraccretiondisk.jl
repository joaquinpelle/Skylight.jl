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

function set_emitter_four_velocity!(vector, position, model::BosonStarAccretionDisk, coord_system::SphericalClass)



end

function get_local_total_flux() = thermal_emission_total_flux()

end


function temperature(positions, )
end
function build_interpolator(temperature)



end

function temperature(position, model::BosonStarAccretionDisk)

    t, r, θ, φ = position





function read_temperature_file(model::BosonStarAccretionDisk)

    temperature = readdlm(model.temperature_file, ' ', Float64, '\n')

end



    
    

