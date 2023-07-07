abstract type AbstractPostProcessCache end

@with_kw mutable struct ImagePlanePostProcessCache{S,M} <: AbstractPostProcessCache
    observer_metric::Matrix{Float64} = zeros(4,4)
    emitter_metric::Matrix{Float64} = zeros(4,4)
    observer_four_velocity::Vector{Float64} = zeros(4)
    emitter_four_velocity::Vector{Float64} = zeros(4)
    spacetime_cache::S
    model_cache::M
end

@with_kw mutable struct PinholeCameraPostProcessCache{S,M} <: AbstractPostProcessCache
    observer_metric::Matrix{Float64} = zeros(4,4)
    emitter_metric::Matrix{Float64} = zeros(4,4)
    observer_four_velocity::Vector{Float64} = zeros(4)
    emitter_four_velocity::Vector{Float64} = zeros(4)
    flux_direction::Vector{Float64} = zeros(4)
    spacetime_cache::S
    model_cache::M
end