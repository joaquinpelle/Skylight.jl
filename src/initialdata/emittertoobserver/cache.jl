function metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations)
    metric_and_tetrad!(cache, position, configurations, opaque_interior_surface_trait(configurations.radiative_model))
    return nothing
end

function metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations, ::IsNotOpaqueInteriorSurface)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    metric, _, time_vector, triad = unpack_views(cache)

    metric!(metric, position, spacetime)
    emitter_four_velocity!(time_vector, position, metric, spacetime, model, coords_top)
    random_triad!(triad, time_vector, metric)
    return nothing
end

function metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations, ::IsOpaqueInteriorSurface)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    metric, metric_inverse, time_vector, triad = unpack_views(cache)

    metric!(metric, position, spacetime)
    metric_inverse!(metric_inverse, position, spacetime)
    emitter_four_velocity!(time_vector, position, metric, spacetime, model, coords_top)
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
    return metric, metric_inverse, time_vector, triad
end