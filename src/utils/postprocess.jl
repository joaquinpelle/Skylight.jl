export get_coordinate_arrays, view_intensities_matrix

function get_coordinate_arrays(configurations)

    image_plane = configurations.image_plane

    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.horizontal_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.horizontal_number_of_nodes

    return range(-0.5*sα, stop=0.5*sα; length=Nα), range(-0.5*sβ, stop=0.5*sβ; length=Nβ)

end

function view_intensities_matrix(output_data, configurations; E_idx)

    image_plane = configurations.image_plane

    NE = length(configurations.observed_energies)
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    return reverse!(reshape(output_data[8+NE+E_idx,:], (Nα, Nβ)), dims=1)

end