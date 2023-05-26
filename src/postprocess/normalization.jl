function normalize_by_image_plane_distance!(array, configurations::AbstractOTEConfigurations)
    image_plane_distance_CGS = geometrized_to_CGS(configurations.image_plane.distance, Dimensions.length, configurations) 
    array ./= image_plane_distance_CGS^2
end

function normalize_by_pixel_area!(array, configurations::AbstractOTEConfigurations)
    pixel_area_CGS = geometrized_to_CGS(pixel_area(configurations.image_plane), Dimensions.area, configurations) 
    array .*= pixel_area_CGS 
end

function rescale_intensities_normalization_at_real_observer!(array, real_observer_distance_CGS, configurations::AbstractOTEConfigurations)
    image_plane_distance_CGS = geometrized_to_CGS(configurations.image_plane.distance, Dimensions.length, configurations) 
    array .*= (image_plane_distance_CGS/real_observer_distance_CGS)^2
end