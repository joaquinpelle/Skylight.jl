function initialize(configurations::VacuumETOConfigurations)
    cache = initial_data_cache(configurations)
    packets = my_zeros(configurations)
    model_cache = allocate_cache(configurations.radiative_model)
    Npp = configurations.number_of_packets_per_point

    for (index, position) in enumerate(initial_positions(configurations))
            packet_index = (index-1)*Npp+1
            @views packets_at_position = packets[:, packet_index:(packet_index+Npp-1)]
            initialize_packets_at_position!(packets_at_position, position, cache, model_cache, configurations)
    end
    return packets
end

function initialize_packets_at_position!(packets_at_position, position, cache, model_cache, configurations)
    
    model = configurations.radiative_model
 
    @views begin
        xμ = packets_at_position[1:4,:]
        kμ = packets_at_position[5:8,:]
    end
    
    metric_and_tetrad!(cache, position, model_cache, configurations)
    
    @views tetrad = cache.tetrad
    
    packets_positions!(xμ, position)
    packets_momenta!(kμ, tetrad, model)
    return nothing
end

function packets_positions!(xμ, position)
    xμ .= position
    return nothing
end

function packets_momenta!(kμ, tetrad, model)
    packets_tetrad_components!(kμ, model)
    coordinate_components_from_tetrad_components!(kμ, tetrad)
    return nothing
end

"""By having Minkowski-null components in the tetrad we guarantee that the momentum is null"""
function packets_tetrad_components!(kμ, model)
    unit_time_components!(kμ)
    trait = opaque_interior_surface_trait(model)
    packets_unit_random_triad_components!(kμ, trait)
    return nothing
end

function packets_unit_random_triad_components!(kμ, ::IsNotOpaqueInteriorSurface)
    @views ki = kμ[2:4,:] 
    random_uniform_points_unit_sphere!(ki, CartesianTopology())
    return nothing
end

"""Sets only positive components along the surface normal"""    
function packets_unit_random_triad_components!(kμ, ::IsOpaqueInteriorSurface)
    @views ki = kμ[2:4,:] 
    random_uniform_points_unit_hemisphere_xaxis!(ki, CartesianTopology())
    return nothing
end