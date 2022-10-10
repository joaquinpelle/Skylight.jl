export initialize_data

function initialize_data(configurations::ETOInitialDataConfigurations)
    
    packets = my_zeros(configurations)
    cache = get_cache(configurations)

    Npp = configurations.number_of_packets_per_point

    index = 1
    for position in get_initial_positions(configurations)

            @views packets_at_position = packets[:, index:(index+Npp-1)]
            initialize_packets_at_position!(packets_at_position, position, cache, configurations)
            index += Npp

    end

    return packets

end

function initialize_packets_at_position!(packets_at_position, position, cache, configurations)

    model = configurations.emission_model

    @views begin
        xμ = packets_at_position[1:4,:]
        kμ = packets_at_position[5:8,:]
    end

    set_metric_and_local_tetrad_in!(cache, position, configurations)
    
    @views tetrad = cache.tetrad
    
    set_packets_position!(xμ, position)
    set_packets_momenta!(kμ, tetrad, model)

end

function set_metric_and_local_tetrad_in!(cache, position, configurations)

    spacetime = configurations.spacetime
    model = configurations.emission_model
    coord_system = coordinate_system_kind(spacetime)
    
    metric, metric_inverse, time_vector, triad = unpack_views(cache)

    set_metric!(metric, position, spacetime)
    set_metric_inverse!(metric_inverse, position, spacetime)
    set_local_four_velocity!(time_vector, position, metric, model, coord_system)     
    set_local_triad!(triad, position, time_vector, metric, metric_inverse, model, coord_system)

end

function set_packets_position!(xμ, position)
    xμ .= position
end

function set_packets_momenta!(kμ, tetrad, model)
    
    set_tetrad_components!(kμ, model)
    set_coordinate_components_from_tetrad_components!(kμ, tetrad)

end

function set_tetrad_components!(kμ, model)
    
    # By having Minkowski-null components in the tetrad we guarantee that the momentum is null

    set_unit_time_component!(kμ)
    set_unit_random_triad_components!(kμ, model)

end

function set_unit_time_component!(kμ)
    kμ[1,:] .= 1.0
end

function set_unit_random_triad_components!(kμ, model)
    
    @views ki = kμ[2:4,:] 
    random_uniform_points_unit_sphere!(ki, CartesianKind())

end

function set_unit_random_triad_components!(kμ, model::OpaqueInteriorSurfaceEmissionModel)
    
    @views ki = kμ[2:4,:] 

    random_uniform_points_unit_hemisphere_xaxis!(ki, CartesianKind())

end