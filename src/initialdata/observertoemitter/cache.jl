#ImagePlaneCache
function unpack_views(cache::ImagePlaneCache)
    @views begin
        metric = cache.metric
        vector = cache.vector
    end
    return metric, vector
end

function dump_∂t_in!(cache::ImagePlaneCache)
    cache.vector = ∂t()
    return nothing
end

#PinholeCameraCache
function unpack_views(cache::PinholeCameraCache)
    @views begin
        metric = cache.metric
        vector = cache.tetrad[:,1]
        triad  = cache.tetrad[:,2:4]
    end
    return metric, vector, triad
end

function set_metric_and_tetrad!(cache::PinholeCameraCache, position, configurations::AbstractConfigurations)
    spacetime = configurations.spacetime
    coords_top = coordinates_topology(spacetime)
    
    metric, time_vector, triad = unpack_views(cache)
    
    set_metric_in!(metric, position, spacetime)
    set_static_four_velocity_in!(time_vector, metric)
    set_spherical_like_triad_in!(triad, position, time_vector, metric, coords_top)
    return nothing
end