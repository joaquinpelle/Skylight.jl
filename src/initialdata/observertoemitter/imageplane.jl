<<<<<<< HEAD:src/initialdata/observertoemitter.jl
function get_initial_data(configurations::AbstractOTEConfigurations)
=======
function get_initial_data(image_plane::ImagePlane, configurations::AbstractOTEConfigurations)
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
    rays = my_zeros(configurations)
    cache = get_initial_data_cache(configurations)

    dump_∂t_in!(cache)
    
    index = 1
    for initial_time in observed_times(configurations) 
        for pixel_coordinates in get_pixel_coordinates(image_plane) 
            @views ray = rays[1:8, index]
            initialize_single!(ray, initial_time, pixel_coordinates, configurations, cache)
            index += 1
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
    space_momentum .= space_momentum_from(pixel_coordinates,image_plane,coords_top)

    dump_metric_in!(cache,position,spacetime)
    set_null_ingoing_past_directed!(momentum,cache)
<<<<<<< HEAD:src/initialdata/observertoemitter.jl
end

function get_space_position_from(pixel_coordinates, image_plane::ImagePlane, ::CartesianTopology)
=======
    return nothing
end

function space_position_from(pixel_coordinates, image_plane::ImagePlane, ::CartesianTopology)
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
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

<<<<<<< HEAD:src/initialdata/observertoemitter.jl
function get_space_position_from(pixel_coordinates, image_plane::ImagePlane, ::SphericalTopology)
=======
function space_position_from(pixel_coordinates, image_plane::ImagePlane, ::SphericalTopology)
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
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

<<<<<<< HEAD:src/initialdata/observertoemitter.jl
function get_space_momentum_from(pixel_coordinates, image_plane::ImagePlane, ::SphericalTopology)
=======
function space_momentum_from(pixel_coordinates, image_plane::ImagePlane, ::SphericalTopology)
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
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

<<<<<<< HEAD:src/initialdata/observertoemitter.jl
function get_space_momentum_from(pixel_coordinates, image_plane::ImagePlane, ::CartesianTopology)
=======
function space_momentum_from(pixel_coordinates, image_plane::ImagePlane, ::CartesianTopology)
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
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
<<<<<<< HEAD:src/initialdata/observertoemitter.jl
=======
    return nothing
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
end

"""returns with unit energy"""
function set_null!(momentum,cache)
    gμν, tμ = unpack_views(cache)
    t2 = norm_squared(tμ,gμν)
    k2 = norm_squared(momentum,gμν)    
    kt = scalar_product(momentum,tμ,gμν)
    α = (-kt + sqrt(kt^2 - k2*t2))/k2
    @. momentum = tμ + α*momentum
<<<<<<< HEAD:src/initialdata/observertoemitter.jl
=======
    return nothing
>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
end

function set_ingoing_past_directed!(momentum)
    @. momentum *= -1
<<<<<<< HEAD:src/initialdata/observertoemitter.jl
end
=======
    return nothing
end



>>>>>>> pinhole:src/initialdata/observertoemitter/imageplane.jl
