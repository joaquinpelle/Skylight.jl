export get_observed_bolometric_fluxes
export get_observed_specific_fluxes
export get_q_factors

function get_observed_bolometric_fluxes(initial_data, output_data, configurations::VacuumOTEConfigurations)
    
    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    nrays = number_of_initial_conditions(configurations)
    observed_bolometric_fluxes = zeros(nrays)

    for i in 1:nrays

        #view ui and uf
        @views begin 

            ui = initial_data[1:8,i]
            uf = output_data[1:8,i]
    
        end

        dump_metrics_and_emitter_four_velocity_in!(cache, ui, uf, configurations)
        
        q = get_q_factor(ui, uf, cache)
                
        emitted_bolometric_flux = get_emitted_bolometric_flux(uf, cache, configurations)
        
        observed_bolometric_fluxes[i] = q^4*emitted_bolometric_flux

    end

    return observed_bolometric_fluxes

end

function get_emitted_bolometric_flux(uf, cache, configurations::VacuumOTEConfigurations)

    spacetime = configurations.spacetime
    model = configurations.model
    coord_system = coordinate_system_class(spacetime)

    @views begin
        
        position = uf[1:4]
        momentum = uf[5:8]

        metric = cache.emitter_metric
        emitter_four_velocity = cache.emitter_four_velocity
    
    end

    #The only difference in the ETO scheme here should be the minus sign in front of the final momentum, right?
    emitted_bolometric_flux = get_emitted_bolometric_flux(position, -momentum, emitter_four_velocity, metric, spacetime, model, coord_system)

    return emitted_bolometric_flux

end

function get_q_factor(ui, uf, cache)

    @views begin

        ki = ui[5:8]
        kf = uf[5:8]
    
    end
    
    q = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)

    return q

end

function get_observed_specific_fluxes(q_factors, output_data, configurations::VacuumOTEConfigurations, frequencies)

    #get_local_emission(rays, model)
    #return Q**3*Fem

end



function get_q_factors(initial_data, output_data, configurations::VacuumOTEConfigurations)

    spacetime = configurations.spacetime
    model = configurations.model
    coord_system = coordinate_system_class(spacetime)
    nrays = number_of_initial_conditions(configurations)

    observer_metric = zeros(4,4)
    emitter_metric = zeros(4,4)
    
    observer_four_velocity = ∂t()
    emitter_four_velocity = zeros(4)

    q_factors = zeros(nrays)

    for i in 1:nrays

        @views begin 

            initial_position = initial_data[1:4,i]
            initial_momentum = initial_data[5:8,i]
            final_position = output_data[1:4,i]
            final_momentum = output_data[5:8,i]
    
        end

        set_metric!(observer_metric, initial_position, spacetime) 
        set_metric!(emitter_metric, final_position, spacetime)

        set_emitter_four_velocity!(emitter_four_velocity, final_position, emitter_metric, spacetime, model, coord_system)

        q_factors[i] = scalar_product(initial_momentum, observer_four_velocity, observer_metric) /
                       scalar_product(final_momentum, emitter_four_velocity, emitter_metric)
        
    end

    return q_factors

end

function get_emitted_fluxes_comoving_frame(output_data, configurations::VacuumOTEConfigurations)