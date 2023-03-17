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

function get_observed_bolometric_fluxes(rays, model::BosonStarAccretionDisk)

    #get_Q_factors
    #get_local_emission
    #return Q**4*Fem

end


function get_observed_specific_fluxes(rays, model::BosonStarAccretionDisk)



end



function temperature(positions, )
function build_interpolator(temperature)



end

function temperature(position, model::BosonStarAccretionDisk)

    t, r, θ, φ = position





function read_temperature_file(model::BosonStarAccretionDisk)

    temperature = readdlm(model.temperature_file, ' ', Float64, '\n')

end



    
    

