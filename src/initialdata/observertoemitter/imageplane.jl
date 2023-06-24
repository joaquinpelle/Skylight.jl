function initialize(image_plane::ImagePlane, configurations::AbstractOTEConfigurations; chunks_per_thread::Int=2)
    rays = my_zeros(configurations)
    Npx = number_of_pixels(image_plane)
    for (it, initial_time) in enumerate(observation_times(configurations)) 
        # Break the work into chunks. More chunks per thread has better load balancing but more overhead
        nchunks = div(Npx, chunks_per_thread*nthreads())
        chunks = Iterators.partition(enumerate(camera_grid(image_plane)), nchunks)
        # Map over the chunks, creating an array of spawned tasks
        map(chunks) do chunk
            Threads.@spawn begin
                cache = initial_data_cache(configurations)
                for (ipx, pixel_coordinates) in chunk 
                    index = (it-1)*Npx+ipx
                    @views ray = rays[1:8, index]
                    initialize_single!(ray, initial_time, pixel_coordinates, configurations, cache)
                end
            end
        end
    end
    return rays
end

function initialize_single!(ray, initial_time, pixel_coordinates, configurations, cache)
    @views begin 
        position = ray[1:4]
        space_position = ray[2:4]
        momentum = ray[5:8]
        space_momentum = ray[6:8]
    end
    
    spacetime = configurations.spacetime
    image_plane = configurations.camera
    coords_top = coordinates_topology(spacetime)
    
    ray[1] = initial_time  
    space_position .= space_position_from(pixel_coordinates,image_plane,coords_top)

    metric!(cache.metric,position,spacetime)
    static_four_velocity!(cache)
    
    space_momentum .= space_momentum_from(pixel_coordinates,image_plane,coords_top)
    set_null_ingoing_past_directed!(momentum,cache)
    return nothing
end

function space_position_from(pixel_coordinates, image_plane::ImagePlane, ::CartesianTopology)
    α,β = pixel_coordinates
    ξ = image_plane.observer_inclination_in_radians
    d = image_plane.distance
    
    sinξ = sin(ξ)
    cosξ = cos(ξ)

    x =-β*cosξ+d*sinξ
    y = α
    z = β*sinξ+d*cosξ

    return [x, y, z]
end

function space_position_from(pixel_coordinates, image_plane::ImagePlane, ::SphericalTopology)
    α,β = pixel_coordinates
    ξ = image_plane.observer_inclination_in_radians
    d = image_plane.distance
    
    sinξ = sin(ξ)
    cosξ = cos(ξ)
    
    r = sqrt(α^2+β^2+d^2)
    θ = acos((d*cosξ+β*sinξ)/r)
    φ = atan(α,(d*sinξ-β*cosξ))

    return [r, θ, φ]
end

function space_momentum_from(pixel_coordinates, image_plane::ImagePlane, ::SphericalTopology)
    α,β = pixel_coordinates
    ξ = image_plane.observer_inclination_in_radians
    d = image_plane.distance
    
    r = sqrt(α^2+β^2+d^2)
    
    sinξ = sin(ξ)
    cosξ = cos(ξ)

    kr =  d/r
    kθ = (-cosξ + d/r^2*(d*cosξ+β*sinξ))/sqrt(r^2-(d*cosξ+β*sinξ)^2)
    kφ =  -α*sinξ/(α^2+(d*sinξ-β*cosξ)^2)  

    return [kr, kθ, kφ]
end

function space_momentum_from(pixel_coordinates, image_plane::ImagePlane, ::CartesianTopology)
    ξ = image_plane.observer_inclination_in_radians
    kx = sin(ξ)
    ky = 0.0
    kz = cos(ξ)
    return [kx, ky, kz]
end

"""the input time component must be zero for this to work """
function set_null_ingoing_past_directed!(momentum,cache)
    set_null!(momentum,cache)
    set_ingoing_past_directed!(momentum)
    return nothing
end

"""returns with unit energy"""
function set_null!(momentum,cache)
    gμν, tμ = unpack_views(cache)
    t2 = norm_squared(tμ,gμν)
    k2 = norm_squared(momentum,gμν)    
    kt = scalar_product(momentum,tμ,gμν)
    α = (-kt + sqrt(kt^2 - k2*t2))/k2
    @. momentum = tμ + α*momentum
    return nothing
end

function set_ingoing_past_directed!(momentum)
    @. momentum *= -1
    return nothing
end