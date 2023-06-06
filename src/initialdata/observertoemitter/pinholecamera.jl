"""
    Sets the bundle of rays at the camera position. The rays are initialized on a tetrad which has a static four-velocity and 
    a spherical-like spatial triad. The rays have unit energy in this tetrad. 
    See docs/pinholecamera.md for more details"""
function initialize(camera::PinholeCamera, configurations::AbstractOTEConfigurations)
    cache = initial_data_cache(camera)
    @views tetrad = cache.tetrad
    
    rays = my_zeros(configurations)    
    position = camera.position
    @views begin
        xμ = rays[1:4,:]
        kμ = rays[5:8,:]
    end
    metric_and_tetrad!(cache, position, configurations)
    rays_position!(xμ, position)
    rays_momenta!(kμ, tetrad, camera)
    return rays
end

function rays_position!(xμ, position)
    xμ .= position
    return nothing
end

function rays_momenta!(kμ, tetrad, camera::PinholeCamera)
    rays_tetrad_components!(kμ, camera)
    coordinate_components_from_tetrad_components!(kμ, tetrad)
    return nothing
end

"""By having Minkowski-null components in the tetrad we guarantee that the momentum is null"""
function rays_tetrad_components!(kμ, camera::PinholeCamera)
    negative_unit_time_components!(kμ)
    rays_triad_components!(kμ, camera)
    return nothing
end

function rays_triad_components!(kμ, camera::PinholeCamera)
    @views space_kμ = kμ[2:4,:]
    i = 1 
    for (α, β) in camera_grid(camera)
        space_kμ[1,i] = cos(α)*cos(β)
        space_kμ[2,i] = sin(α)*cos(β)
        space_kμ[3,i] = sin(β)
        i += 1
    end
    return nothing        
end