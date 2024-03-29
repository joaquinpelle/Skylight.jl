#ImagePlaneCache
function static_four_velocity!(cache::ImagePlaneCache)
    static_four_velocity!(cache.vector, cache.metric)
    return nothing
end

function unpack_views(cache::ImagePlaneCache)
    @views begin
        metric = cache.metric
        vector = cache.vector
        scache = cache.spacetime_cache
    end
    return metric, vector, scache
end

#PinholeCameraCache

function tetrad!(cache::PinholeCameraCache, camera, spacetime::AbstractSpacetime)
    metric!(cache.metric, camera.position, spacetime, cache.spacetime_cache)
    v = default_four_velocity(camera, spacetime)
    tetrad!(cache, camera.position, v, spacetime)
    return nothing
end

function tetrad!(cache::PinholeCameraCache,
    position,
    four_velocity,
    spacetime::AbstractSpacetime)
    coords_top = coordinates_topology(spacetime)
    metric, time_vector, triad, _ = unpack_views(cache)
    time_vector .= four_velocity
    spherical_like_triad!(triad, position, time_vector, metric, coords_top)
    return nothing
end

function unpack_views(cache::PinholeCameraCache)
    @views begin
        metric = cache.metric
        vector = cache.tetrad[:, 1]
        triad = cache.tetrad[:, 2:4]
        scache = cache.spacetime_cache
    end
    return metric, vector, triad, scache
end
