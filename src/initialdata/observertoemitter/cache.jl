#ImagePlaneCache
function static_four_velocity!(cache::ImagePlaneCache)
    static_four_velocity!(cache.vector, cache.metric)
    return nothing
end

function unpack_views(cache::ImagePlaneCache)
    @views begin
        metric = cache.metric
        vector = cache.vector
    end
    return metric, vector
end

#PinholeCameraCache

function tetrad!(cache::PinholeCameraCache, position, spacetime::AbstractSpacetime)
    coords_top = coordinates_topology(spacetime)
    metric, time_vector, triad = unpack_views(cache)
    static_four_velocity!(time_vector, metric)
    spherical_like_triad!(triad, position, time_vector, metric, coords_top)
    return nothing
end

function unpack_views(cache::PinholeCameraCache)
    @views begin
        metric = cache.metric
        vector = cache.tetrad[:,1]
        triad  = cache.tetrad[:,2:4]
    end
    return metric, vector, triad
end