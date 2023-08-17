"""
    Sets the bundle of rays at the camera position. The rays are initialized on a tetrad which has a static four-velocity and 
    a spherical-like spatial triad. The rays have unit energy in this tetrad. 
    See docs/pinholecamera.md for more details"""
function initialize(camera::PinholeCamera, configurations::AbstractOTEConfigurations)
    spacetime = configurations.spacetime
    position = camera.position
    four_velocity!(camera, spacetime)
    cache = initial_data_cache(configurations)
    metric!(cache.metric, position, spacetime,cache.spacetime_cache)
    has_lorentzian_signature(cache.metric) || throw(ArgumentError("The metric signature is not Lorentzian."))
    tetrad!(cache, position, camera.four_velocity, spacetime)
    rays = my_zeros(configurations)    
    @views begin
        xμ = rays[1:4,:]
        kμ = rays[5:8,:]
        tetrad = cache.tetrad
    end
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
    @inbounds begin
        pixels = (collect∘grid)(camera)
        @threads for i in eachindex(pixels)
            α, β = pixels[i]
            space_kμ[1,i] = cos(α)*cos(β)
            space_kμ[2,i] = sin(α)*cos(β)
            space_kμ[3,i] = sin(β)
        end
    end
    return nothing        
end