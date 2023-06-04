function set_metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations)
    set_metric_and_tetrad!(cache, position, configurations, opaque_interior_surface_trait(configurations.radiative_model))
end

function set_metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations, ::IsNotOpaqueInteriorSurface)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    metric, time_vector, triad = unpack_views(cache)

    set_metric!(metric, position, spacetime)
    set_emitter_four_velocity!(time_vector, position, metric, spacetime, model, coords_top)
    set_random_triad!(triad, time_vector, metric)
end

function set_metric_and_tetrad!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations, ::IsOpaqueInteriorSurface)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    metric, metric_inverse, time_vector, triad = unpack_views(cache)

    set_metric!(metric, position, spacetime)
    set_metric_inverse!(metric_inverse, position, spacetime)
    set_emitter_four_velocity!(time_vector, position, metric, spacetime, model, coords_top)
    set_surface_adapted_triad!(triad, time_vector, metric, metric_inverse, position, model, coords_top)
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