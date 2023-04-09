export get_observed_bolometric_fluxes
export get_observed_specific_fluxes

function get_observed_bolometric_fluxes(initial_data, output_data, configurations::VacuumOTEConfigurations)
    
    spacetime = configurations.spacetime
    model = configurations.model
    coord_system = coordinate_system_class(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    observed_bolometric_fluxes = zeros(Nrays)

    for i in 1:Nrays

        @views begin 

            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]

        end

        dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coord_system)
        
        q = get_energies_quotient(ki, kf, cache)

        #Check: the only difference with the ETO scheme here should be the minus sign in front of the final momentum
        emitted_bolometric_flux = get_emitted_bolometric_flux(pf, -kf, cache.emitter_four_velocity, cache.metric, spacetime, model, coord_system)

        observed_bolometric_fluxes[i] = q^4*emitted_bolometric_flux

    end

    return observed_bolometric_fluxes

end

function get_observed_specific_fluxes(initial_data, output_data, observed_energies, configurations::VacuumOTEConfigurations)

    spacetime = configurations.spacetime
    model = configurations.model
    coord_system = coordinate_system_class(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    NE = length(observed_energies)

    observed_specific_fluxes = zeros(NE, Nrays)

    for i in 1:Nrays

        @views begin 
        
            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]

        end

        dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coord_system)
        
        q = get_energies_quotient(ki, kf, cache)

        for j in 1:NE

            emitted_energy = observed_energies[j]/q

            #Check: the only difference with the ETO scheme here should be the minus sign in front of the final momentum
            emitted_specific_flux = get_emitted_specific_flux(pf, -kf, emitted_energy, cache.emitter_four_velocity, cache.metric, spacetime, model, coord_system)
            observed_specific_fluxes[j, i] = q^3*emitted_specific_flux
        
        end

    end
    
    return observed_specific_fluxes

end


function get_energies_quotient(ki, kf, cache)
    
    q = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)

    return q

end

