@with_kw mutable struct ImagePlaneCache 
    metric::Matrix{Float64} = zeros(4,4)
    vector::Vector{Float64} = zeros(4)
end

@with_kw mutable struct PinholeCameraCache 
    metric::Matrix{Float64} = zeros(4,4)
    tetrad::Vector{Float64} = zeros(4,4)
end

@with_kw mutable struct ETOInitialDataCache
    metric::Matrix{Float64} = zeros(4,4)
    metric_inverse::Matrix{Float64} = zeros(4,4)
    tetrad::Matrix{Float64} = zeros(4,4)
end