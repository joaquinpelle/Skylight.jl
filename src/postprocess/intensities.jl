function observed_bolometric_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations)
    #Returns intensities in CGS and energy quotiens

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    q = zeros(Nrays)
    Iobs = zeros(Nrays)

    @inbounds begin
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

            dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coords_top)
            q[i] = energies_quotient(ki, kf, cache)
            Iobs[i] = q[i]^4*emitted_bolometric_intensity(pf, -kf, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
        end
    end

    normalize_by_image_plane_distance!(Iobs, configurations)
    return Iobs, q
end

function observed_specific_intensities(initial_data, output_data, observed_energies_CGS, configurations::VacuumOTEConfigurations)
    #Returns intensities in CGS
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    NE = length(observed_energies)

    q = zeros(Nrays)
    Iobs = zeros(NE, Nrays)
    @inbounds begin
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

            dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coords_top)
            q[i] = energies_quotient(ki, kf, cache)
            
            for j in 1:NE
                emitted_energy = observed_energies_CGS[j]/q[i]
                Iobs[j, i] = q[i]^3*emitted_specific_intensity(pf, -kf, emitted_energy, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
            end
        end
    end

    normalize_by_image_plane_distance!(Iobs, configurations)
    return Iobs, q
end