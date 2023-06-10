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

function observer_four_velocity!(cache, observer_four_velocity) 
    if observer_four_velocity===nothing
        static_four_velocity!(cache.observer_four_velocity, cache.observer_metric)
    else
        is_timelike(observer_four_velocity, cache.observer_metric) || throw(ArgumentError("The observer four-velocity is not timelike."))
        cache.observer_four_velocity .= observer_four_velocity
    end
end

function flux_direction!(cache, flux_direction, camera, spacetime)
    if flux_direction===nothing
        cache.flux_direction .= default_flux_direction(camera, spacetime)
    else
        is_spacelike(flux_direction, cache.observer_metric) || throw(ArgumentError("The flux direction is not spacelike."))
        cache.flux_direction = flux_direction
    end
end

function emitter_metric_and_four_velocity!(cache::PinholeCameraPostProcessCache, final_position, spacetime, model, coords_top)
    metric!(cache.emitter_metric, final_position, spacetime)
    emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coords_top)
    return nothing
end
