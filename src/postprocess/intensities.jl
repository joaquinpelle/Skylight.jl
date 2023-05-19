function get_observed_bolometric_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations)
    
    #Returns intensities in CGS

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coord_system = coordinate_system_class(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    observed_bolometric_intensities = zeros(Nrays)

    for i in 1:Nrays

        @views begin 

            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]

        end

        #The difference with the ETO scheme here should be the minus sign in front of the final momentum
        #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end

        dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coord_system)
        
        q = get_energies_quotient(ki, kf, cache)

        emitted_bolometric_intensity = get_emitted_bolometric_intensity(pf, -kf, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coord_system)

        observed_bolometric_intensities[i] = q^4*emitted_bolometric_intensity
        
    end

    set_intensities_normalization_at_image_plane!(observed_bolometric_intensities, configurations)

    return observed_bolometric_intensities

end

function get_observed_specific_intensities(initial_data, output_data, observed_energies_CGS, configurations::VacuumOTEConfigurations)

    #Returns intensities in CGS

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coord_system = coordinate_system_class(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    NE = length(observed_energies)

    observed_specific_intensities = zeros(NE, Nrays)

    for i in 1:Nrays

        @views begin 
        
            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]

        end

        #The difference with the ETO scheme here should be the minus sign in front of the final momentum
        #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end

        dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coord_system)
        
        q = get_energies_quotient(ki, kf, cache)

        for j in 1:NE

            emitted_energy = observed_energies_CGS[j]/q

            emitted_specific_intensity = get_emitted_specific_intensity(pf, -kf, emitted_energy, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coord_system)
            observed_specific_intensities[j, i] = q^3*emitted_specific_intensity
        
        end

    end

    set_intensities_normalization_at_image_plane!(observed_specific_intensities, configurations)
    
    return observed_specific_intensities

end

function get_energies_quotient(ki, kf, cache)
    
    q = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)

    return q

end

function set_intensities_normalization_at_image_plane!(intensities, configurations)

    image_plane_distance_CGS = geometrized_to_CGS(configurations.image_plane.distance, Dimensions.length, configurations) 
    
    intensities ./= image_plane_distance_CGS^2

end

function rescale_intensities_normalization_at_real_observer!(intensities, real_observer_distance_CGS, configurations)

    image_plane_distance_CGS = geometrized_to_CGS(configurations.image_plane.distance, Dimensions.length, configurations) 
    
    intensities .*= (image_plane_distance_CGS/real_observer_distance_CGS)^2

end