#AbstractPostProcessCache

#ImagePlanePostProcessCache
function metrics_and_four_velocities!(cache::ImagePlanePostProcessCache, initial_position, final_position, spacetime, model, coords_top)
    metric!(cache.observer_metric, initial_position, spacetime) 
    metric!(cache.emitter_metric, final_position, spacetime)
    static_four_velocity!(cache.observer_four_velocity, cache.observer_metric)
    emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coords_top)
    return nothing
end

function unpack_views(cache::ImagePlanePostProcessCache)
    @views begin
        observer_metric = cache.observer_metric
        emitter_metric = cache.emitter_metric
        observer_four_velocity = cache.observer_four_velocity
        emitter_four_velocity = cache.emitter_four_velocity
    end
    return observer_metric, emitter_metric, observer_four_velocity, emitter_four_velocity
end

#PinholeCameraPostProcessCache
function observer_metric!(cache::PinholeCameraPostProcessCache, position, spacetime) 
    metric!(cache.observer_metric, position, spacetime)
    return nothing
end

function emitter_metric_and_four_velocity!(cache::PinholeCameraPostProcessCache, final_position, spacetime, model, coords_top)
    metric!(cache.emitter_metric, final_position, spacetime)
    emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coords_top)
    return nothing
end