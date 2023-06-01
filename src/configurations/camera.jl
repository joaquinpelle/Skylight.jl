#AbstractCamera methods
function get_pixel_coordinates(camera::AbstractCamera)
    sα, sβ = sides(camera)
    Nα, Nβ = numbers_of_nodes_per_side(camera)

    horizontal_coordinates = range(-0.5*sα, stop=0.5*sα; length=Nα)
    vertical_coordinates = range(-0.5*sβ,0.5*sβ; length=Nβ)

    return Iterators.product(horizontal_coordinates,vertical_coordinates)
end

function grid_spacing(camera::AbstractCamera)
    sα, sβ = sides(camera)
    Nα, Nβ = numbers_of_nodes_per_side(camera)

    dα = sα/(Nα-1)
    dβ = sβ/(Nβ-1)
    return dα, dβ
end

function number_of_nodes(camera::AbstractCamera)
    Nα, Nβ = numbers_of_nodes_per_side(camera)
    return Nα*Nβ
end

function numbers_of_nodes_per_side(camera::AbstractCamera)
    Nα = camera.horizontal_number_of_nodes
    Nβ = camera.vertical_number_of_nodes
    return Nα, Nβ
end

#PinholeCamera methods
function sides(camera::PinholeCamera)
    sα = camera.horizontal_aperture
    sβ = camera.vertical_aperture
    return sα, sβ
end

time(camera::PinholeCamera) = camera.position[1]

#ImagePlane methods
function area(image_plane::ImagePlane)
    Nα, Nβ = numbers_of_nodes_per_side(image_plane)
    dA = pixel_area(image_plane)
    return Nα*Nβ*dA
end

function pixel_area(image_plane::ImagePlane)
    dα, dβ = grid_spacing(image_plane)
    return dα*dβ 
end

function sides(image_plane::ImagePlane)
    sα = image_plane.horizontal_side
    sβ = image_plane.vertical_side
    return sα, sβ
end

function get_rmax(image_plane::ImagePlane) 
    d = image_plane.distance
    sα, sβ = sides(image_plane)
    return 1.1*sqrt(d^2 + sα^2 + sβ^2)
end

time(image_plane::ImagePlane) = 0.0