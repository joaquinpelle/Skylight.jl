spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; Emin::Number, Emax::Number, NE::Integer, kwargs...) = spectrum(initial_data, output_data, configurations, range(Emin, stop=Emax, length=NE); kwargs...)
spectrum(initial_data, output_data, configurations::VacuumOTEConfigurationsm, energies; kwargs...) = spectrum(initial_data, output_data, configurations, configurations.camera, energies; kwargs...)

function spectrum(initial_data, 
                output_data, 
                configurations::VacuumOTEConfigurations,
                camera::ImagePlane,
                energies; 
                kwargs...)
    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    Iobs, _ = observed_specific_intensities(initial_data, output_data, configurations, energies; kwargs...)
    Fobs = fluxes(Iobs, camera)
    return sum(Fobs, dims=2)
end

function spectrum(initial_data, 
                output_data, 
                configurations::VacuumOTEConfigurations,
                camera::PinholeCamera,
                energies; 
                kwargs...)
    
    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    Iobs, _ = observed_specific_intensities(initial_data, output_data, configurations, energies; kwargs...)
    Fobs = fluxes(Iobs, camera, initial_data, output_data, configurations.spacetime; kwargs...)
    return sum(Fobs, dims=2)
end