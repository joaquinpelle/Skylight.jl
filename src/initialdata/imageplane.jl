export ImagePlane

@with_kw struct ImagePlane

    observer_distance :: Float64
    observer_inclination_in_degrees :: Float64
    horizontal_side_image_plane :: Float64
    vertical_side_image_plane :: Float64
    horizontal_number_of_nodes :: Int64
    vertical_number_of_nodes :: Int64
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)

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

