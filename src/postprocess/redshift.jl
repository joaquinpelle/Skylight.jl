function energies_quotients(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations;
    kwargs...)
    energies_quotients(initial_data,
        output_data,
        configurations,
        configurations.camera;
        kwargs...)
end

function energies_quotients(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    camera::ImagePlane;
    tasks_per_thread::Int = 2)
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    nrays = size(initial_data, 2)
    q = zeros(nrays)
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = div(nrays, nthreads() * tasks_per_thread)
    chunks = Iterators.partition(1:nrays, chunk_size)
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    @sync map(chunks) do chunk
        Threads.@spawn begin
            cache = postprocess_cache(configurations)
            for i in chunk
                @views begin
                    pi = initial_data[1:4, i]
                    ki = initial_data[5:8, i]

                    pf = output_data[1:4, i]
                    kf = output_data[5:8, i]
                end
                #The difference with the ETO method here should be the minus sign in front of the final momentum
                #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
                if !is_final_position_at_source(pf, spacetime, model)
                    continue
                end
                metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
                q[i] = energies_quotient(ki, kf, cache)
            end
        end
    end
    return q
end

function energies_quotients(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    camera::PinholeCamera;
    observer_four_velocity = nothing,
    tasks_per_thread::Int = 2)
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    nrays = size(initial_data, 2)
    q = zeros(nrays)
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = div(nrays, nthreads() * tasks_per_thread)
    chunks = Iterators.partition(1:nrays, chunk_size)
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    @sync map(chunks) do chunk
        Threads.@spawn begin
            cache = postprocess_cache(configurations)
            observer_metric!(cache, camera.position, spacetime)
            observer_four_velocity!(cache, observer_four_velocity)
            for i in chunk
                @views begin
                    ki = initial_data[5:8, i]
                    pf = output_data[1:4, i]
                    kf = output_data[5:8, i]
                end
                #The difference with the ETO method here should be the minus sign in front of the final momentum
                #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
                if !is_final_position_at_source(pf, spacetime, model)
                    continue
                end
                emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
                q[i] = energies_quotient(ki, kf, cache)
            end
        end
    end
    return q
end

function energies_quotient(ki, kf, cache)
    return scalar_product(ki, cache.observer_four_velocity, cache.observer_metric) /
           scalar_product(kf, cache.rest_frame_four_velocity, cache.emitter_metric)
end

"""
Assuming all photons are emitted with unit initial energy
"""
function energies_quotients(data::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractAccretionDisk)
    coords_top = coordinates_topology(spacetime)
    nrays = size(data, 2)
    q = zeros(nrays)
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = div(nrays, Threads.nthreads()*2)
    chunks = Iterators.partition(1:nrays, chunk_size)
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    @sync map(chunks) do chunk
        Threads.@spawn begin
            g = zeros(4,4)
            u = zeros(4)
            spacetime_cache = allocate_cache(spacetime)
            model_cache = allocate_cache(model)
            for i in chunk
                @views begin 
                    position = data[1:4,i]
                    momentum = data[5:8,i]
                end
                metric!(g, position, spacetime, spacetime_cache)
                rest_frame_four_velocity!(u, position, g, spacetime, model, coords_top, spacetime_cache, model_cache)
                q[i] = -Skylight.scalar_product(u,momentum,g)
            end
        end
    end
    return q
end