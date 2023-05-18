function my_zeros(configurations::NonVacuumConfigurations)

    NE = length(configurations.observed_energies)
    
    return zeros(8+2*NE, number_of_initial_conditions(configurations))

end

my_zeros(configurations::VacuumConfigurations) = zeros(8, number_of_initial_conditions(configurations))

get_observed_times(configurations::OTEConfigurations) = configurations.observed_times

get_initial_data_cache(::OTEConfigurations) = OTEInitialDataCache()
get_initial_data_cache(::ETOConfigurations) = ETOInitialDataCache()

get_postprocess_cache(::OTEConfigurations) = OTEPostProcessCache()

function get_initial_positions(configurations::ETOConfigurations)
    
    times = zero_times(configurations)
    space_positions = get_space_positions(configurations)
    
    return eachcol([times'; space_positions])

end

function get_space_positions(configurations::ETOConfigurations)
    
    npoints = configurations.number_of_points

    coord_system = coordinate_system_class(configurations.spacetime)
    space_positions = get_space_positions(npoints, configurations.radiative_model, coord_system)

    return space_positions

end

function zero_times(configurations::ETOConfigurations)
    
    npoints = configurations.number_of_points
    return repeat([0.0],npoints)

end

function number_of_initial_conditions(configurations::OTEConfigurations)
     
    number_of_times = length(configurations.observed_times)
    
    return number_of_nodes(configurations.image_plane)*number_of_times 
    
end

function number_of_initial_conditions(configurations::ETOConfigurations)
    
    number_of_points = configurations.number_of_points
    number_of_packets_per_point = configurations.number_of_packets_per_point

    return number_of_points*number_of_packets_per_point
    
end

function get_pixel_coordinates(image_plane::ImagePlane)
    
    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.vertical_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes

    horizontal_coordinates = range(-0.5*sα, stop=0.5*sα; length=Nα)
    vertical_coordinates = range(-0.5*sβ,0.5*sβ; length=Nβ)

    return Iterators.product(horizontal_coordinates,vertical_coordinates)

end

function number_of_nodes(image_plane::ImagePlane)
    return image_plane.horizontal_number_of_nodes*image_plane.vertical_number_of_nodes
end

function pixel_area(image_plane::ImagePlane)
    
    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.vertical_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    dα = sα/(Nα-1)
    dβ = sβ/(Nβ-1)
    
    return dα*dβ 

end

function area(image_plane::ImagePlane)

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes

    dA = pixel_area(image_plane)

    return Nα*Nβ*dA

end