#AbstractCamera methods
function get_pixel_coordinates(camera::AbstractCamera)
    α, β = get_pixel_coordinates_vectors(camera)
    return Iterators.product(α, β)
end

function get_pixel_coordinates_vectors(camera::AbstractCamera)
    sα, sβ = sides(camera)
    Nα, Nβ = numbers_of_pixels_per_side(camera)
    dα, dβ = grid_spacing(camera)
    return range(-0.5*sα+0.5*dα, stop=0.5*sα-0.5*dα; length=Nα), range(-0.5*sβ+0.5*dβ, stop=0.5*sβ-0.5*dβ; length=Nβ)
end

function grid_spacing(camera::AbstractCamera)
    sα, sβ = sides(camera)
    Nα, Nβ = numbers_of_pixels_per_side(camera)
    dα = sα/Nα
    dβ = sβ/Nβ
    return dα, dβ
end

function number_of_pixels(camera::AbstractCamera)
    Nα, Nβ = numbers_of_pixels_per_side(camera)
    return Nα*Nβ
end

function numbers_of_pixels_per_side(camera::AbstractCamera)
    Nα = camera.horizontal_number_of_nodes
    Nβ = camera.vertical_number_of_nodes
    return Nα, Nβ
end

#PinholeCamera methods
function solid_angle(camera::PinholeCamera)
    sα, sβ = sides(camera)
    return 2*sα*sin(sβ/2)
end

function sides(camera::PinholeCamera)
    sα = camera.horizontal_aperture_in_radians
    sβ = camera.vertical_aperture_in_radians
    return sα, sβ
end

time(camera::PinholeCamera) = camera.position[1]

get_initial_data_cache(::PinholeCamera) = PinholeCameraCache()

#ImagePlane methods
function area(image_plane::ImagePlane)
    sα, sβ = sides(image_plane)
    return sα*sβ
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

get_initial_data_cache(::ImagePlane) = ImagePlaneCache()