export ImagePlane
export OTEInitialDataConfigurations
export initialize_OTE

@with_kw struct ImagePlane

    observer_distance :: Float64
    observer_inclination_in_degrees :: Float64
    horizontal_side_image_plane :: Float64
    vertical_side_image_plane :: Float64
    horizontal_number_of_nodes :: Int32
    vertical_number_of_nodes :: Int32
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)

end

@with_kw struct OTEInitialDataConfigurations{T<:Spacetime} 
    
    spacetime::T
    image_plane::ImagePlane 
    initial_times::Vector{Float64}

end

function initialize_OTE(configurations::OTEInitialDataConfigurations)

    rays = zero_rays_on_grid(configurations)

    container = zeros(4,5)
    dump_∂t_in!(container)
    
    index = 1
    
    for initial_time in get_initial_times(configurations)
        for pixel_coordinates in get_pixel_coordinates(configurations)

            @views ray = rays[index,1:8]
            initialize_single!(ray, initial_time, pixel_coordinates, configurations, container)
            index += 1

        end
    end

    return rays

end

get_initial_times(configurations) = configurations.initial_times

function get_pixel_coordinates(configs::OTEInitialDataConfigurations)
    
    image_plane = configs.image_plane
    
    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.vertical_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes

    horizontal_coordinates = range(-0.5*sα, stop=0.5*sα; length=Nα)
    vertical_coordinates = range(-0.5*sβ,0.5*sβ; length=Nβ)

    return Iterators.product(horizontal_coordinates,vertical_coordinates)

end

function initialize_single!(ray, initial_time, pixel_coordinates, configurations, container)
    
    @views begin 
        space_position = ray[2:4]
        momentum = ray[5:8]
        space_momentum = ray[6:8]
        position = ray[1:4]
    end
    
    spacetime = configurations.spacetime
    coord_system = spacetime.asymptotic_coordinate_system
    image_plane = configurations.image_plane


    ray[1] = initial_time
    space_position .= get_space_position_from(pixel_coordinates,image_plane,coord_system)
    space_momentum .= get_space_momentum_from(pixel_coordinates,image_plane,coord_system)

    dump_metric_in!(container,position,spacetime)
    set_null_ingoing_past_directed!(momentum,container)

end

function get_space_position_from(pixel_coordinates, image_plane::ImagePlane, coord_system::CartesianCoordinates)

    α,β = pixel_coordinates
    ξ = image_plane.observer_inclination_in_radians
    d = image_plane.observer_distance
    
    sinξ = sin(ξ)
    cosξ = cos(ξ)

    x =-β*cosξ+d*sinξ
    y = α
    z = β*sinξ+d*cosξ

    return [x, y, z]

end

function get_space_position_from(pixel_coordinates, image_plane::ImagePlane, coord_system::SphericalCoordinates)

    α,β = pixel_coordinates
    ξ = image_plane.observer_inclination_in_radians
    d = image_plane.observer_distance
    
    sinξ = sin(ξ)
    cosξ = cos(ξ)
    
    r = sqrt(α^2+β^2+d^2)
    θ = acos((d*cosξ+β*sinξ)/r)
    φ = atan(α,(d*sinξ-β*cosξ))

    return [r, θ, φ]

end

function get_space_momentum_from(pixel_coordinates, image_plane::ImagePlane, coord_system::SphericalCoordinates)

    α,β = pixel_coordinates
    ξ = image_plane.observer_inclination_in_radians
    d = image_plane.observer_distance
    
    r = sqrt(α^2+β^2+d^2)
    
    sinξ = sin(ξ)
    cosξ = cos(ξ)

    kr =  d/r
    kθ = (-cosξ + d/r^2*(d*cosξ+β*sinξ))/sqrt(r^2-(d*cosξ+β*sinξ)^2)
    kφ =  -α*sinξ/(α^2+(d*sinξ-β*cosξ)^2)  

    return [kr, kθ, kφ]

end

function get_space_momentum_from(pixel_coordinates, image_plane::ImagePlane, coord_system::CartesianCoordinates)
    
    ξ = image_plane.observer_inclination_in_radians
    
    kx = sin(ξ)
    ky = 0.0
    kz = cos(ξ)

    return [kx, ky, kz]

end

function dump_∂t_in!(container)
    container[:,5] = ∂t()
end

function dump_metric_in!(container,position,spacetime::Spacetime)

    pars = spacetime.parameters
    metric! = spacetime.metric!
    
    @views g = container[:,1:4]
    metric!(g,position,pars)
    
end

function set_null_ingoing_past_directed!(momentum,container)
    
    """the input time component must be zero for this to work """

    set_null!(momentum,container)
    set_ingoing_past_directed!(momentum)

end

function set_null!(momentum,container)
    
    """returns with unit energy"""

    @views begin
        gμν = container[:,1:4]
        tμ = container[:,5]
    end

    t2 = norm_squared(tμ,gμν)
    k2 = norm_squared(momentum,gμν)    
    kt = scalar_product(momentum,tμ,gμν)

    α = (-kt + sqrt(kt^2 - k2*t2))/k2

    @. momentum = tμ + α*momentum
    
end

function set_ingoing_past_directed!(momentum)
    @. momentum *= -1
end

function pixel_area(image_plane::ImagePlane)
    
    sα = image_plane.horizontal_side_image_plane
    sβ = image_plane.vertical_side_image_plane
    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    
    dα = sα/(Nα-1)
    dβ = sβ/(Nβ-1)
    
    return dα*dβ 

end

function area(image_plane::ImagePlane)

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes

    dA = pixel_area(image_plane)

    return Nα*Nβ*dA

end

function zero_rays_on_grid(configurations::OTEInitialDataConfigurations)
    
    image_plane = configurations.image_plane

    Nα = image_plane.horizontal_number_of_nodes
    Nβ = image_plane.vertical_number_of_nodes
    Nt = length(configurations.initial_times)

    return zeros(Nα*Nβ*Nt,8)
    
end