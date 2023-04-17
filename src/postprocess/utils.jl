export get_coordinate_arrays, view_intensities_grid

function get_coordinate_arrays(configurations::OTEConfigurations)

    image_plane = configurations.image_plane

    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.horizontal_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.horizontal_number_of_nodes

    return range(-0.5*sα, stop=0.5*sα; length=Nα), range(-0.5*sβ, stop=0.5*sβ; length=Nβ)

end

function view_intensities_grid(output_data, configurations::NonVacuumConfigurations, E_idx)

    image_plane = configurations.image_plane

    NE = length(configurations.observed_energies)
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    @views intensities_grid = reverse!(reshape(output_data[8+NE+E_idx,:], (Nα, Nβ)), dims=1)
    
    return intensities_grid

end

function view_intensities_grid(intensities, configurations::VacuumOTEConfigurations)

    image_plane = configurations.image_plane

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    @views intensities_grid = reverse!(reshape(intensities, (Nα, Nβ)), dims=1)
    
    return intensities_grid

end

function view_intensities_grid(intensities, configurations::VacuumOTEConfigurations, E_idx)

    image_plane = configurations.image_plane

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    @views intensities_grid = reverse!(reshape(intensities[E_idx,:], (Nα, Nβ)), dims=1)
    
    return intensities_grid

end