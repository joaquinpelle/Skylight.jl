#AbstractPostProcessCache

#ImagePlanePostProcessCache
function set_metrics_and_four_velocities!(cache::ImagePlanePostProcessCache, initial_position, final_position, spacetime, model, coords_top)
    set_metric!(cache.observer_metric, initial_position, spacetime) 
    set_metric!(cache.emitter_metric, final_position, spacetime)
    set_static_four_velocity!(cache.observer_four_velocity, cache.observer_metric)
    set_emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coords_top)
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
function set_observer_metric!(cache::PinholeCameraPostProcessCache, position, spacetime) 
    set_metric!(cache.observer_metric, position, spacetime)
    return nothing
end

function set_observer_four_velocity!(cache::PinholeCameraPostProcessCache, observer_four_velocity) 
    normalize_timelike!(observer_four_velocity, cache.observer_metric)
    cache.observer_four_velocity .= observer_four_velocity
    return nothing
end

function set_observer_four_velocity!(cache::PinholeCameraPostProcessCache, ::Nothing)
    set_static_four_velocity!(cache.observer_four_velocity, cache.observer_metric)
    return nothing
end

function set_surface_normal!(cache::PinholeCameraPostProcessCache, configurations, normal)
    normalize_spacelike!(normal, cache.observer_metric)
    cache.surface_element_normal .= normal
    return nothing
end

function set_surface_normal!(cache::PinholeCameraPostProcessCache, configurations, ::Nothing)
    cache.surface_element_normal .= default_normal(configurations.camera, configurations)
    return nothing
end

function set_emitter_metric_and_four_velocity!(cache::PinholeCameraPostProcessCache, final_position, spacetime, model, coords_top)
    set_metric!(cache.emitter_metric, final_position, spacetime)
    set_emitter_four_velocity!(cache.emitter_four_velocity, final_position, cache.emitter_metric, spacetime, model, coords_top)
    return nothing
end