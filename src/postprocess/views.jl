@kwdispatch grid_view

@kwmethod function grid_view(output_data, configurations::NonVacuumOTEConfigurations; energy_index)
    NE = length(configurations.observation_energies)
    1 <= energy_index <= NE || throw(ArgumentError("energy_index must be between 1 and $NE"))
    Nα, Nβ = numbers_of_pixels_per_side(configurations.camera)
    @views grid = reshape(output_data[8+NE+energy_index,:], (Nα, Nβ))
    return grid
end

function grid_view(vector, configurations::VacuumOTEConfigurations)
    Nα, Nβ = numbers_of_pixels_per_side(configurations.camera)
    @views grid = reshape(vector, (Nα, Nβ))
    return grid
end

@kwmethod function grid_view(intensities, configurations::VacuumOTEConfigurations; energy_index)
    Nα, Nβ = numbers_of_pixels_per_side(configurations.camera)
    @views intensities_grid = reshape(intensities[energy_index,:], (Nα, Nβ))
    return intensities_grid
end