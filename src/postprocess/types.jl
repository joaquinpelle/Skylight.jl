abstract type AbstractPostProcessCache end

@with_kw mutable struct OTEPostProcessCache <: AbstractPostProcessCache
    
    observer_metric::Matrix{Float64} = zeros(4,4)
    emitter_metric::Matrix{Float64} = zeros(4,4)
    observer_four_velocity::Vector{Float64} = zeros(4)
    emitter_four_velocity::Vector{Float64} = zeros(4)

end