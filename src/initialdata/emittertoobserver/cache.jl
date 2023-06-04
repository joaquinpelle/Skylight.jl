function unpack_views(cache::ETOInitialDataCache)
    @views begin
        metric = cache.metric
        metric_inverse = cache.metric_inverse
        vector = cache.tetrad[:,1]
        triad  = cache.tetrad[:,2:4]
    end
    return metric, metric_inverse, vector, triad
end

function dump_metric_and_tetrad_in!(cache::ETOInitialDataCache, position, configurations::AbstractConfigurations)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    
    dump_metric_in!(cache, position, spacetime)
    dump_metric_inverse_in!(cache, position, spacetime)
    dump_tetrad_in!(cache, position, spacetime, model, coords_top)
end

function dump_metric_inverse_in!(cache::ETOInitialDataCache, position, spacetime::AbstractSpacetime)
    set_metric_inverse!(cache.metric_inverse, position, spacetime)
end

function dump_tetrad_in!(cache::ETOInitialDataCache, position, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology)
    dump_emitter_four_velocity_in!(cache, position, spacetime, model, coords_top)     
    dump_triad_in!(cache, position, model, coords_top)
end

function dump_emitter_four_velocity_in!(cache::ETOInitialDataCache, position, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology)
    @views time_vector = cache.tetrad[:,1]
    set_emitter_four_velocity!(time_vector, position, cache.metric, spacetime, model, coords_top)
end

function dump_triad_in!(cache::ETOInitialDataCache, position, ::AbstractRadiativeModel, ::AbstractCoordinatesTopology)
    metric, metric_inverse, time_vector, triad = unpack_views(cache)
    set_triad!(triad, time_vector, metric)
end

function dump_triad_in!(cache::ETOInitialDataCache, position, model::AbstractSurfaceEmissionModel, coords_top::AbstractCoordinatesTopology)
    metric, metric_inverse, time_vector, triad = unpack_views(cache)
    set_surface_adapted_triad!(triad, time_vector, position, metric, metric_inverse, model, coords_top)
end