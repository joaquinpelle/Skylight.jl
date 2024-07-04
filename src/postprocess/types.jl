abstract type AbstractPostProcessCache end

@with_kw mutable struct ImagePlanePostProcessCache{S, M} <: AbstractPostProcessCache
    observer_metric::Matrix{Float64} = zeros(4, 4)
    emitter_metric::Matrix{Float64} = zeros(4, 4)
    observer_four_velocity::Vector{Float64} = zeros(4)
    rest_frame_four_velocity::Vector{Float64} = zeros(4)
    spacetime_cache::S
    model_cache::M
end

@with_kw mutable struct PinholeCameraPostProcessCache{S, M} <: AbstractPostProcessCache
    observer_metric::Matrix{Float64} = zeros(4, 4)
    emitter_metric::Matrix{Float64} = zeros(4, 4)
    observer_four_velocity::Vector{Float64} = zeros(4)
    rest_frame_four_velocity::Vector{Float64} = zeros(4)
    flux_direction::Vector{Float64} = zeros(4)
    spacetime_cache::S
    model_cache::M
end

@with_kw mutable struct ETOPostProcessCache{S, M} <: AbstractPostProcessCache
    emitter_metric::Matrix{Float64} = zeros(4, 4)
    emitter_metric_inverse::Matrix{Float64} = zeros(4, 4)
    rest_frame_four_velocity::Vector{Float64} = zeros(4)
    surface_normal::Vector{Float64} = zeros(4)
    spacetime_cache::S
    model_cache::M
end