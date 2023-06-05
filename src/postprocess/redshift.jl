function energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations, camera::ImagePlane)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    cache = postprocess_cache(camera)
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
        metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
    end
    return q
end

function energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations, camera::PinholeCamera; observer_four_velocity=nothing)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)

    if observer_four_velocity===nothing
        static_four_velocity!(cache.observer_four_velocity, cache.observer_metric)
    else
        @assert is_timelike(observer_four_velocity, cache.observer_metric) "The observer four-velocity is not timelike."
        cache.observer_four_velocity .= observer_four_velocity
    end
    
    Nrays = number_of_initial_conditions(configurations)
    q = zeros(Nrays)
    for i in 1:Nrays
        @views begin 
            ki = initial_data[5:8,i]
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]
        end
        #The difference with the ETO scheme here should be the minus sign in front of the final momentum
        #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end
        emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
    end
    return q
end

function energies_quotient(ki, kf, cache)
    return scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)
end