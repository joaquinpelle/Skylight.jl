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
    set_packets_tetrad_components!(kμ, model)
    set_coordinate_components_from_tetrad_components!(kμ, tetrad)
end

"""By having Minkowski-null components in the tetrad we guarantee that the momentum is null"""
function set_packets_tetrad_components!(kμ, model)
    set_unit_time_components!(kμ)
    trait = opaque_interior_surface_trait(model)
    set_packets_unit_random_triad_components!(kμ, trait)
end

function set_packets_unit_random_triad_components!(kμ, ::IsNotOpaqueInteriorSurface)
    @views ki = kμ[2:4,:] 
    random_uniform_points_unit_sphere!(ki, CartesianTopology())
end

"""Sets only positive components along the surface normal"""    
function set_packets_unit_random_triad_components!(kμ, ::IsOpaqueInteriorSurface)
    @views ki = kμ[2:4,:] 
    random_uniform_points_unit_hemisphere_xaxis!(ki, CartesianTopology())
end