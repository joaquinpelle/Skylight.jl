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

function dump_∂t_in!(cache)
    cache.vector = ∂t()
end

function dump_metric_in!(cache, position, spacetime)
    set_metric!(cache.metric, position, spacetime)
end

function dump_metric_inverse_in!(cache, position, spacetime)
    set_metric_inverse!(cache.metric_inverse, position, spacetime)
end

function dump_metric_and_tetrad_in!(cache, position, configurations)

    spacetime = configurations.spacetime
    model = configurations.emission_model
    coord_system = coordinate_system_kind(spacetime)
    
    dump_metric_in!(cache, position, spacetime)
    dump_metric_inverse_in!(cache, position, spacetime)
    dump_tetrad_in!(cache, position, model, coord_system)

end

function dump_tetrad_in!(cache, position, model, coord_system)
    
    dump_model_four_velocity_in!(cache, position, model, coord_system)     
    dump_triad_in!(cache, position, model, coord_system)

end

function dump_model_four_velocity_in!(cache, position, model, coord_system)
    @views time_vector = cache.tetrad[:,1]
    set_model_four_velocity!(time_vector, position, cache.metric, model, coord_system)
end


function dump_triad_in!(cache, position, model, coord_system)

    metric, metric_inverse, time_vector, triad = unpack_views(cache)
    set_triad!(triad, time_vector, metric)

end

function dump_triad_in!(cache, position, model::SurfaceEmissionModel, coord_system)

    metric, metric_inverse, time_vector, triad = unpack_views(cache)
    set_surface_adapted_triad!(triad, time_vector, position, metric, metric_inverse, model, coord_system)

end
