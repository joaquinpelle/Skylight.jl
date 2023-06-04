"""
    Builds a bundle of rays set on a tetrad at the camera position. The tetrad has the normalized ∂t as time vector and the central direction
    projected orthogonally to the time vector as the first spatial vector. The rays have unit energy in this tetrad. 
    See docs/pinholecamera.md for more details"""
function get_initial_data(camera::PinholeCamera, configurations::AbstractOTEConfigurations)
    cache = get_initial_data_cache(camera)
    @views tetrad = cache.tetrad
    
    rays = my_zeros(configurations)    
    position = copy(camera.position)

    Npixels = number_of_pixels(camera) 
    i = 1
    for t in observation_times(configurations)
        position[1] += t
        @views begin
            xμ = rays[1:4,i:(i-1+Npixels)]
            kμ = rays[5:8,i:(i-1+Npixels)]
        end
        set_metric_and_tetrad!(cache, position, configurations)
        set_rays_position!(xμ, position)
        set_rays_momenta!(kμ, tetrad, camera)
        i += Npixels
    end
    return rays
end

function set_rays_position!(xμ, position)
    xμ .= position
    return nothing
end

function set_rays_momenta!(kμ, tetrad, camera::PinholeCamera)
    set_rays_tetrad_components!(kμ, camera)
    set_coordinate_components_from_tetrad_components!(kμ, tetrad)
    return nothing
end

"""By having Minkowski-null components in the tetrad we guarantee that the momentum is null"""
function set_rays_tetrad_components!(kμ, camera::PinholeCamera)
    set_negative_unit_time_component!(kμ)
    set_rays_triad_components!(kμ, camera)
    return nothing
end

function set_rays_triad_components!(kμ, camera::PinholeCamera)
    @views space_kμ = kμ[2:4,:]
    i = 1 
    for (α, β) in get_pixel_coordinates(camera)
        @. begin
            space_kμ[1,i] = cos(α)*cos(β)
            space_kμ[2,i] = sin(α)*cos(β)
            space_kμ[3,i] = sin(β)
        end
        i += 1
    end
    return nothing        
end