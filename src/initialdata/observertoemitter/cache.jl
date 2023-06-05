#ImagePlaneCache
function metric_and_four_velocity!(cache::ImagePlaneCache, position, spacetime::AbstractSpacetime)
    metric, time_vector = unpack_views(cache)
    metric!(metric, position, spacetime)
    static_four_velocity!(time_vector, metric)
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
function metric_and_tetrad!(cache::PinholeCameraCache, position, configurations::AbstractConfigurations)
    spacetime = configurations.spacetime
    coords_top = coordinates_topology(spacetime)
    metric, time_vector, triad = unpack_views(cache)
    metric!(metric, position, spacetime)
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
