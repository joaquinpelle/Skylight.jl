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
    Nα = camera.horizontal_number_of_pixels
    Nβ = camera.vertical_number_of_pixels
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

function all_pixel_solid_angles(camera::PinholeCamera)
    Nα, _ = numbers_of_pixels_per_side(camera)
    _, vecβ = get_pixel_coordinates_vectors(camera)
    dα, dβ = grid_spacing(camera)

    solid_angles = my_zeros(camera)
    i=1
    for β in vecβ
        solid_angle[i:(i-1+Nα)] .=  2*cos(β)*sin(dβ/2)*dα
        i+=Nα
    end
    return solid_angles
end

function default_tetrad(camera::PinholeCamera, configurations::AbstractOTEConfigurations, t)
    cache = get_initial_data_cache(camera)
    position = coyp(camera.position)
    position[1] += t
    set_metric_and_tetrad!(cache, position, configurations)
    return cache.tetrad
end

default_four_velocity(camera::PinholeCamera, configurations::AbstractOTEConfigurations, t) = default_tetrad(camera, configurations, t)[:,1]
default_normal(camera::PinholeCamera, configurations::AbstractOTEConfigurations, t) = default_tetrad(camera, configurations, t)[:,2]

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