function spectrum(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations;
    Emin::Number,
    Emax::Number,
    NE::Integer,
    kwargs...)
    spectrum(initial_data,
        output_data,
        configurations,
        range(Emin, stop = Emax, length = NE);
        kwargs...)
end
function spectrum(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    energies;
    kwargs...)
    spectrum(initial_data,
        output_data,
        configurations,
        configurations.camera,
        energies;
        kwargs...)
end
function spectrum(initial_data,
    output_data,
    configurations::NonVacuumOTEConfigurations;
    kwargs...)
    spectrum(initial_data, output_data, configurations, configurations.camera; kwargs...)
end

function spectrum(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    camera::ImagePlane,
    energies)
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))
    Fobs, _ = observed_specific_intensities(initial_data,
        output_data,
        configurations,
        energies)
    fluxes!(Fobs, configurations, camera)
    return sum(Fobs, dims = 2)[:, 1]
end

function spectrum(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    camera::PinholeCamera,
    energies;
    observer_four_velocity = nothing,
    flux_direction = nothing,
    tasks_per_thread::Int = 2)
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))
    eight_components(initial_data, output_data) ||
        throw(DimensionMismatch("The initial and output data must have eight components."))
    number_of_pixels(camera) == size(initial_data, 2) ||
        throw(DimensionMismatch("The number of pixels in the camera must be the same as the number of data rows."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    nrays = size(initial_data, 2)
    NE = length(energies)

    Fobs = zeros(NE, nrays)
    dΩ = pixel_solid_angles(camera)
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = div(nrays, nthreads() * tasks_per_thread)
    chunks = Iterators.partition(1:nrays, chunk_size)
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    @sync map(chunks) do chunk
        Threads.@spawn begin
            cache = postprocess_cache(configurations)
            observer_metric!(cache, camera.position, spacetime)
            observer_four_velocity!(cache, observer_four_velocity)
            flux_direction!(cache, flux_direction, camera, spacetime)
            for i in chunk
                @views begin
                    ki = initial_data[5:8, i]
                    pf = output_data[1:4, i]
                    kf = output_data[5:8, i]
                end
                if !is_final_position_at_source(pf, spacetime, model)
                    continue
                end
                emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
                q = energies_quotient(ki, kf, cache)
                for j in axes(energies, 1)
                    emitted_energy = energies[j] / q
                    #The difference with the ETO scheme here should be the minus sign in front of the final momentum
                    #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
                    Fem = rest_frame_specific_intensity(pf,
                        -kf,
                        emitted_energy,
                        cache.rest_frame_four_velocity,
                        cache.emitter_metric,
                        spacetime,
                        model,
                        coords_top,
                        cache.model_cache)
                    nu = scalar_product(ki,
                        cache.observer_four_velocity,
                        cache.observer_metric)
                    nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
                    Fobs[j, i] = q^3 * Fem * nu * nn * dΩ[i]
                end
            end
        end
    end
    return sum(Fobs, dims = 2)[:, 1]
end

function spectrum(initial_data,
    output_data,
    configurations::NonVacuumOTEConfigurations,
    camera::ImagePlane)
    Fobs = observed_specific_intensities(initial_data, output_data, configurations)
    fluxes!(Fobs, configurations, camera)
    return sum(Fobs, dims = 2)[:, 1]
end

function spectrum(initial_data,
    output_data,
    configurations::NonVacuumOTEConfigurations,
    camera::PinholeCamera;
    observer_four_velocity = nothing,
    flux_direction = nothing,
    tasks_per_thread::Int = 2)
    number_of_pixels(camera) == size(initial_data, 2) ||
        throw(DimensionMismatch("The number of pixels in the camera must be the same as the number of data rows."))

    spacetime = configurations.spacetime
    observation_energies = configurations.observation_energies
    nrays = size(initial_data, 2)
    NE = length(observation_energies)
    Fobs = zeros(NE, nrays)
    dΩ = pixel_solid_angles(camera)
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = div(nrays, nthreads() * tasks_per_thread)
    chunks = Iterators.partition(1:nrays, chunk_size)
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    @sync map(chunks) do chunk
        Threads.@spawn begin
            cache = postprocess_cache(configurations)
            observer_metric!(cache, camera.position, spacetime)
            observer_four_velocity!(cache, observer_four_velocity)
            flux_direction!(cache, flux_direction, camera, spacetime)
            for i in chunk
                @views begin
                    ki = initial_data[5:8, i]
                    pf = output_data[1:4, i]
                    kf = output_data[5:8, i]
                end
                observer_rest_frame_energy = scalar_product(ki,
                    cache.observer_four_velocity,
                    cache.observer_metric)
                nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
                nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
                @. Fobs[:, i] = (observation_energies * observer_rest_frame_energy)^3 *
                                output_data[(9 + NE):end, i]
                Fobs[:, i] *= nu * nn * dΩ[i]
            end
        end
    end
    return sum(Fobs, dims = 2)[:, 1]
end
