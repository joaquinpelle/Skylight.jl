abstract type InitialDataCache end

@with_kw mutable struct OTEInitialDataCache <: InitialDataCache
    metric::Matrix{Float64} = zeros(4,4)
    vector::Vector{Float64} = zeros(4)
end

@with_kw mutable struct ETOInitialDataCache <: InitialDataCache
    metric::Matrix{Float64} = zeros(4,4)
    metric_inverse::Matrix{Float64} = zeros(4,4)
    tetrad::Matrix{Float64} = zeros(4,4)
end

function dump_vector_in!(cache, vector)
    cache.vector = vector
end

function dump_∂t_in!(cache)
    cache.vector = ∂t()
end

function dump_metric_in!(cache, position, spacetime::Spacetime)
    set_metric!(cache.metric, position, spacetime)
end

function dump_metric_inverse_in!(cache, position, spacetime::Spacetime)
    set_metric_inverse!(cache.metric_inverse, position, spacetime)
end

function unpack_views(cache::OTEInitialDataCache)

    @views begin
        metric = cache.metric
        vector = cache.vector
    end

    return metric, vector

end

function unpack_views(cache::ETOInitialDataCache)

    @views begin
        metric = cache.metric
        metric_inverse = cache.metric_inverse
        vector = cache.tetrad[:,1]
        triad  = cache.tetrad[:,2:4]
    end

    return metric, metric_inverse, vector, triad

end