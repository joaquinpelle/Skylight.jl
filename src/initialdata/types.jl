abstract type AbstractInitialDataCache end

@with_kw mutable struct ImagePlaneCache <: AbstractInitialDataCache 
    metric::Matrix{Float64} = zeros(4,4)
    vector::Vector{Float64} = zeros(4)
end

@with_kw mutable struct PinholeCameraCache <: AbstractInitialDataCache 
    metric::Matrix{Float64} = zeros(4,4)
    tetrad::Matrix{Float64} = zeros(4,4)
end

@with_kw mutable struct ETOInitialDataCache <: AbstractInitialDataCache
    metric::Matrix{Float64} = zeros(4,4)
    metric_inverse::Matrix{Float64} = zeros(4,4)
    tetrad::Matrix{Float64} = zeros(4,4)
end