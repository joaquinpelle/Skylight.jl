#AbstractCamera methods
function grid(camera::AbstractCamera)
    α, β = axes_ranges(camera)
    return Iterators.product(α, β)
end

function axes_ranges(camera::AbstractCamera)
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

number_of_initial_conditions(camera::PinholeCamera) = number_of_pixels(camera)  

function solid_angle(camera::PinholeCamera)
    sα, sβ = sides(camera)
    return 2*sα*sin(sβ/2)
end

function sides(camera::PinholeCamera)
    sα = camera.horizontal_aperture_in_radians
    sβ = camera.vertical_aperture_in_radians
    return sα, sβ
end

function pixel_solid_angles(camera::PinholeCamera)
    Nα, Nβ = numbers_of_pixels_per_side(camera)
    _, vecβ = axes_ranges(camera)
    dα, dβ = grid_spacing(camera)

    solid_angles = zeros(Nα*Nβ)
    for (i,β) in enumerate(vecβ)
        index = (i-1)*Nα+1
        solid_angles[index:(index-1+Nα)] .=  2*cos(β)*sin(dβ/2)*dα
    end
    return solid_angles
end

function four_velocity!(camera::PinholeCamera, spacetime::AbstractSpacetime)
    if all(camera.four_velocity .== 0.0)
        camera.four_velocity .= static_four_velocity(metric(camera.position, spacetime))
    end 
end

function default_tetrad(camera::PinholeCamera, spacetime::AbstractSpacetime)
    cache = PinholeCameraCache(spacetime)
    metric!(cache.metric, camera.position, spacetime, cache.spacetime_cache)
    tetrad!(cache, camera.position, camera.four_velocity, spacetime)
    return cache.tetrad
end

function default_four_velocity(camera::PinholeCamera, spacetime::AbstractSpacetime)
    g = zeros(4,4)
    spacetime_cache = allocate_cache(spacetime)
    metric!(g, camera.position, spacetime, spacetime_cache)
    return static_four_velocity(g)
end

default_flux_direction(camera::PinholeCamera, spacetime::AbstractSpacetime) = default_tetrad(camera, spacetime)[:,2]
max_radius(camera, spacetime) = 1.1*radius(camera.position, spacetime)

initial_data_cache(::PinholeCamera, spacetime::AbstractSpacetime) = PinholeCameraCache(spacetime)
PinholeCameraCache(spacetime::AbstractSpacetime) = PinholeCameraCache(spacetime_cache = allocate_cache(spacetime))

function postprocess_cache(::PinholeCamera, spacetime::AbstractSpacetime, model::AbstractRadiativeModel) 
    return PinholeCameraPostProcessCache(spacetime_cache = allocate_cache(spacetime), 
                                        model_cache = allocate_cache(model))
end

#ImagePlane methods

number_of_initial_conditions(camera::ImagePlane) = number_of_pixels(camera)*number_of_times(camera)  

function number_of_times(camera::ImagePlane)
    return length(camera.observation_times)
end

function area(camera::ImagePlane)
    sα, sβ = sides(camera)
    return sα*sβ
end

function pixel_area(camera::ImagePlane)
    dα, dβ = grid_spacing(camera)
    return dα*dβ 
end

function sides(camera::ImagePlane)
    sα = camera.horizontal_side
    sβ = camera.vertical_side
    return sα, sβ
end

function max_radius(camera::ImagePlane, ::AbstractSpacetime) 
    d = camera.distance
    sα, sβ = sides(camera)
    return 1.1*sqrt(d^2 + sα^2 + sβ^2)
end

initial_data_cache(::ImagePlane, spacetime::AbstractSpacetime) = ImagePlaneCache(spacetime)
ImagePlaneCache(spacetime::AbstractSpacetime) = ImagePlaneCache(spacetime_cache = allocate_cache(spacetime))

function postprocess_cache(::ImagePlane, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)
    return ImagePlanePostProcessCache(spacetime_cache = allocate_cache(spacetime), 
                               model_cache = allocate_cache(model))
end