function set_metrics_and_four_velocities!(cache::OTEPostProcessCache, initial_position, final_position, spacetime, model, coords_top)
    set_metric!(cache.observer_metric, initial_position, spacetime) 
    set_metric!(cache.emitter_metric, final_position, spacetime)
    set_static_four_velocity!(cache.observer_four_velocity, cache.observer_metric)
    set_emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coords_top)
end

function unpack_views(cache::OTEPostProcessCache)
    @views begin
        observer_metric = cache.observer_metric
        emitter_metric = cache.emitter_metric
        observer_four_velocity = cache.observer_four_velocity
        emitter_four_velocity = cache.emitter_four_velocity
    end
    return observer_metric, emitter_metric, observer_four_velocity, emitter_four_velocity
end