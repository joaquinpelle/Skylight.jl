function ETOInitialDataCache(spacetime::AbstractSpacetime, model::AbstractRadiativeModel)
    scache = allocate_cache(spacetime)
    mcache = allocate_cache(model)
    return ETOInitialDataCache(spacetime_cache=scache, model_cache=mcache)
end

function metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations)
    metric_and_tetrad!(cache, position, configurations, opaque_interior_surface_trait(configurations.radiative_model))
    return nothing
end

function metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations, ::IsNotOpaqueInteriorSurface)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    metric, _, time_vector, triad, scache, mcache = unpack_views(cache)

    metric!(metric, position, spacetime, scache)
    emitter_four_velocity!(time_vector, position, metric, spacetime, model, coords_top, mcache)
    random_triad!(triad, time_vector, metric)
    return nothing
end

function metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations, ::IsOpaqueInteriorSurface)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    metric, metric_inverse, time_vector, triad, scache, mcache = unpack_views(cache)

    metric!(metric, position, spacetime, scache)
    metric_inverse!(metric_inverse, position, spacetime)
    emitter_four_velocity!(time_vector, position, metric, spacetime, model, coords_top, mcache)
    surface_adapted_triad!(triad, time_vector, metric, metric_inverse, position, model, coords_top)
    return nothing
end

function unpack_views(cache::ETOInitialDataCache)
    @views begin
        metric = cache.metric
        metric_inverse = cache.metric_inverse
        time_vector = cache.tetrad[:,1]
        triad  = cache.tetrad[:,2:4]
    end
    return metric, metric_inverse, time_vector, triad, cache.spacetime_cache, cache.model_cache
end