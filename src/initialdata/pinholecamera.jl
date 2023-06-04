
"""
    Builds a bundle of rays set on a tetrad at the camera position. The tetrad has the normalized ∂t as time vector and the central direction
    projected orthogonally to the time vector as the first spatial vector. The rays have unit energy in this tetrad"""
function get_initial_data(camera::PinholeCamera, configurations::AbstractOTEConfigurations)
    
    cache = get_initial_data_cache(camera)
    rays = my_zeros(configurations)

    @views begin
        xμ = rays[1:4,:]
        kμ = rays[5:8,:]
    end

    position = camera.position
    dump_metric_and_tetrad_in!(cache, position, configurations)
    
    @views tetrad = cache.tetrad

    set_rays_position!(xμ, position)
    set_rays_momenta!(kμ, tetrad, camera)

    index = 1
    for initial_time in observed_times(configurations) #TODO in the case of the PinholeCamera this should add to the position[1] of camera
        for pixel_coordinates in get_pixel_coordinates(camera) 
            @views ray = rays[1:8, index]
            initialize_single!(ray, initial_time, pixel_coordinates, spacetime, camera, cache)
            index += 1
        end
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
    α, β = get_pixel_coordinates_vectors(camera)
    @. begin
        space_kμ[1,:] = cos(α)*cos(β)
        space_kμ[2,:] =-sin(α)*cos(β)
        space_kμ[3,:] = sin(β)
    end
    return nothing        
end
