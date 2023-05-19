function get_initial_data(configurations::VacuumETOConfigurations)
    
    cache = get_initial_data_cache(configurations)
    packets = my_zeros(configurations)

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

    model = configurations.radiative_model

    @views begin
        xμ = packets_at_position[1:4,:]
        kμ = packets_at_position[5:8,:]
    end

    dump_metric_and_tetrad_in!(cache, position, configurations)
    
    @views tetrad = cache.tetrad
    
    set_packets_position!(xμ, position)
    set_packets_momenta!(kμ, tetrad, model)

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
    random_uniform_points_unit_sphere!(ki, CartesianClass())

end

function set_unit_random_triad_components!(kμ, ::OpaqueInteriorSurfaceEmissionModel)
    
    #Sets only positive components along the surface normal    
    
    @views ki = kμ[2:4,:] 
    random_uniform_points_unit_hemisphere_xaxis!(ki, CartesianClass())

end