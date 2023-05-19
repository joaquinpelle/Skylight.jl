function get_pixel_coordinates_vectors(configurations::AbstractOTEConfigurations)

    image_plane = configurations.image_plane

    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.horizontal_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.horizontal_number_of_nodes

    return range(-0.5*sα, stop=0.5*sα; length=Nα), range(-0.5*sβ, stop=0.5*sβ; length=Nβ)

end

@kwdispatch view_intensities_grid

@kwmethod function view_intensities_grid(output_data, configurations::NonVacuumOTEConfigurations; E_idx)

    image_plane = configurations.image_plane

    NE = length(configurations.observed_energies)
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    @views intensities_grid = reshape(output_data[8+NE+E_idx,:], (Nα, Nβ))
    
    return intensities_grid

end

function view_intensities_grid(intensities, configurations::VacuumOTEConfigurations)

    image_plane = configurations.image_plane

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    @views intensities_grid = reshape(intensities, (Nα, Nβ))
    
    return intensities_grid

end

@kwmethod function view_intensities_grid(intensities, configurations::VacuumOTEConfigurations; E_idx)

    image_plane = configurations.image_plane

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    @views intensities_grid = reshape(intensities[E_idx,:], (Nα, Nβ))
    
    return intensities_grid

end