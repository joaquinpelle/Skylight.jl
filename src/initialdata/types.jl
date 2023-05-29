@with_kw struct ImagePlane

    distance::Float64
    observer_inclination_in_degrees::Float64
    horizontal_side::Float64
    vertical_side::Float64
    horizontal_number_of_nodes::Int
    vertical_number_of_nodes::Int
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)

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