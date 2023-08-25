abstract type AbstractInitialDataCache end

@with_kw mutable struct ImagePlaneCache{SC} <: AbstractInitialDataCache
    metric::Matrix{Float64} = zeros(4, 4)
    vector::Vector{Float64} = zeros(4)
    spacetime_cache::SC
end

@with_kw mutable struct PinholeCameraCache{SC} <: AbstractInitialDataCache
    metric::Matrix{Float64} = zeros(4, 4)
    tetrad::Matrix{Float64} = zeros(4, 4)
    spacetime_cache::SC
end

@with_kw mutable struct ETOInitialDataCache{SC, MC} <: AbstractInitialDataCache
    metric::Matrix{Float64} = zeros(4, 4)
    metric_inverse::Matrix{Float64} = zeros(4, 4)
    tetrad::Matrix{Float64} = zeros(4, 4)
    spacetime_cache::SC
    model_cache::MC
end
