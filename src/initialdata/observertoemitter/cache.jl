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

function dump_metric_and_tetrad_in!(cache::PinholeCameraCache, position, configurations::AbstractConfigurations)
    spacetime = configurations.spacetime
    coords_top = coordinates_topology(spacetime)
    dump_metric_in!(cache, position, spacetime)
    dump_tetrad_in!(cache, position, cache.metric, coords_top)
end

function dump_tetrad_in!(cache::PinholeCameraCache, position, metric, coords_top::AbstractCoordinatesTopology)
    dump_static_four_velocity_in!(cache, metric)
    dump_ingoing_central_direction_in!(cache, position, cache.tetrad[:,1], metric, coords_top)
    dump_otrhonormal_dyad_in!(cache, cache.tetrad[:,1], cache.tetrad[:,2], metric)
end

function dump_static_four_velocity_in!(cache::PinholeCameraCache, metric)
    @views time_vector = cache.tetrad[:,1]
    time_vector = ∂t()
    normalize_timelike!(time_vector, metric)
end

function dump_ingoing_central_direction_in!(cache::PinholeCameraCache, position, time_vector, metric, ::CartesianTopology)
    @views space_vector = cache.tetrad[:,2]    
    space_vector[1] = 0.0
    space_vector[2] = -position[2]
    space_vector[3] = -position[3]
    space_vector[4] = -position[4]
    space_vector .= orthogonal_projection(space_vector, time_vector, metric)
    normalize_spacelike!(space_vector, metric)
end

function dump_ingoing_central_direction_in!(cache::PinholeCameraCache, position, time_vector, metric, ::SphericalTopology)
    @views space_vector = cache.tetrad[:,2]
    fill!(space_vector, 0.0)
    space_vector[2] = -1.0
    space_vector .= orthogonal_projection(space_vector, time_vector, metric)
    normalize_spacelike!(space_vector, metric)
end

function dump_otrhonormal_dyad_in!(cache::PinholeCamera, time_vector, space_vector, metric)
    @views begin
        dyad = cache.tetrad[:,3:4]
        dyad_time_components = dyad[1,:]
        dyad_space_components = dyad[2:4,:]
    end
    fill!(dyad_time_components,0.0)
    rand!(dyad_space_components)
    orthonormalize!(dyad, time_vector, space_vector, metric)
end