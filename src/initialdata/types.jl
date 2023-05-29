abstract type AbstractCamera end

@with_kw struct ImagePlane <: AbstractCamera
    distance::Float64
    observer_inclination_in_degrees::Float64
    horizontal_side::Float64
    vertical_side::Float64
    horizontal_number_of_nodes::Int
    vertical_number_of_nodes::Int
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)
end

@with_kw struct PinholeCamera <: AbstractCamera 
    position::Vector{Float64}
    four_velocity::Vector{Float64}
    horizontal_aperture::Float64 #This is distance*cos(horizontal_aperture_angle)
    vertical_aperture::Float64   #This is distance*cos(horizontal_aperture_angle)
    horizontal_number_of_nodes::Int
    vertical_number_of_nodes::Int
    direction::Vector{Float64} = to_center(position, four_velocity)
end

@with_kw mutable struct OTEInitialDataCache
    metric::Matrix{Float64} = zeros(4,4)
    vector::Vector{Float64} = zeros(4)
end

@with_kw mutable struct ETOInitialDataCache
    metric::Matrix{Float64} = zeros(4,4)
    metric_inverse::Matrix{Float64} = zeros(4,4)
    tetrad::Matrix{Float64} = zeros(4,4)
end