function energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations)

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)

    Nrays = number_of_initial_conditions(configurations)
    q = zeros(Nrays)

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

        set_metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
    end
    return q
end

function energies_quotient(ki, kf, cache)
    return scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)
end