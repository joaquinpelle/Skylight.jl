energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations; kwargs...) = energies_quotients(initial_data, output_data, configurations, configurations.camera; kwargs...)

function energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations, camera::ImagePlane)
    
    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    threads_cache = postprocess_cache(camera)
    model_cache = allocate_cache(model)
    Nrays = size(initial_data, 2)
    q = zeros(Nrays)
    @inbounds begin
        @threads for i in axes(initial_data, 2)
            cache = threads_cache[threadid()]
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
            metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top, model_cache)
            q[i] = energies_quotient(ki, kf, cache)
        end
    end
    return q
end

function energies_quotients(initial_data, output_data, configurations::VacuumOTEConfigurations, camera::PinholeCamera; observer_four_velocity=nothing)

    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    threads_cache = postprocess_cache(camera)
    model_cache = allocate_cache(model)
    for cache in threads_cache
        observer_metric!(cache, camera.position, spacetime)
        observer_four_velocity!(cache, observer_four_velocity)
    end

    Nrays = size(initial_data, 2)
    q = zeros(Nrays)
    @inbounds begin
        @threads for i in axes(initial_data, 2)
            cache = threads_cache[threadid()]
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
            emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top, model_cache)
            q[i] = energies_quotient(ki, kf, cache)
        end
    end
    return q
end

function energies_quotient(ki, kf, cache)
    return scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
        scalar_product(kf, cache.emitter_four_velocity, cache.emitter_metric)
end