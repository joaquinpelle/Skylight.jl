function spectrum(initial_data, 
                output_data, 
                configurations::VacuumOTEConfigurations,
                camera::ImagePlane; 
                Emin, 
                Emax, 
                NE,
                kwargs...)
    
    energies = range(Emin, stop=Emax, length=NE)
    Iobs, _ = observed_specific_intensities(initial_data, output_data, configurations, energies; kwargs...)
    Fobs = fluxes(Iobs, camera)
    return sum(Fobs, dims=2)
end

function spectrum(initial_data, 
                output_data, 
                configurations::VacuumOTEConfigurations,
                camera::PinholeCamera; 
                Emin, 
                Emax, 
                NE,
                kwargs...)
    
    energies = range(Emin, stop=Emax, length=NE)
    Iobs, _ = observed_specific_intensities(initial_data, output_data, configurations, energies; kwargs...)
    Fobs = fluxes(Iobs, camera, initial_data, output_data, configurations.spacetime; kwargs...)
    return sum(Fobs, dims=2)
end