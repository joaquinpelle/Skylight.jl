@kwdispatch view_as_grid

@kwmethod function view_as_grid(output_data, configurations::NonVacuumOTEConfigurations; E_idx)
    NE = length(configurations.observation_energies)
    Nα, Nβ = numbers_of_pixels_per_side(configurations.camera)
    @views grid = reshape(output_data[8+NE+E_idx,:], (Nα, Nβ))
    return grid
end

function view_as_grid(vector, configurations::VacuumOTEConfigurations)
    Nα, Nβ = numbers_of_pixels_per_side(configurations.camera)
    @views grid = reshape(vector, (Nα, Nβ))
    return grid
end

@kwmethod function view_as_grid(intensities, configurations::VacuumOTEConfigurations; E_idx)
    Nα, Nβ = numbers_of_pixels_per_side(configurations.camera)
    @views intensities_grid = reshape(intensities[E_idx,:], (Nα, Nβ))
    return intensities_grid
end