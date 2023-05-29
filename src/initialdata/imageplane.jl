function get_pixel_coordinates(image_plane::ImagePlane)
    sα, sβ = sides(image_plane)
    Nα, Nβ = numbers_of_nodes_per_side(image_plane)

    horizontal_coordinates = range(-0.5*sα, stop=0.5*sα; length=Nα)
    vertical_coordinates = range(-0.5*sβ,0.5*sβ; length=Nβ)

    return Iterators.product(horizontal_coordinates,vertical_coordinates)
end

function area(image_plane::ImagePlane)
    Nα, Nβ = numbers_of_nodes_per_side(image_plane)
    dA = pixel_area(image_plane)
    return Nα*Nβ*dA
end

function pixel_area(image_plane::ImagePlane)
    dα, dβ = grid_spacing(image_plane)
    return dα*dβ 
end

function grid_spacing(image_plane::ImagePlane)
    sα, sβ = sides(image_plane)
    Nα, Nβ = numbers_of_nodes_per_side(image_plane)

    dα = sα/(Nα-1)
    dβ = sβ/(Nβ-1)
    return dα, dβ
end

function number_of_nodes(image_plane::ImagePlane)
    Nα, Nβ = numbers_of_nodes_per_side(image_plane)
    return Nα*Nβ
end

function numbers_of_nodes_per_side(image_plane::ImagePlane)
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    return Nα, Nβ
end

function sides(image_plane::ImagePlane)
    sα = image_plane.horizontal_side
    sβ = image_plane.vertical_side
    return sα, sβ
end
