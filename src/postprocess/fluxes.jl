function fluxes!(F::AbstractArray, configurations, ::ImagePlane)
    normalize_by_pixel_area!(F, configurations)
    normalize_by_camera_distance!(F, configurations)
end

function fluxes!(F::AbstractArray,
    configurations,
    camera::PinholeCamera,
    initial_data::AbstractMatrix,
    output_data::AbstractMatrix;
    observer_four_velocity = nothing,
    flux_direction = nothing,
    tasks_per_thread::Int = 2)
    fluxes_checks(F, initial_data, output_data, configurations)
    dΩ = pixel_solid_angles(camera)
    itr = axes(F, ndims(F))
    tmap(fluxes_task,
        itr,
        F,
        initial_data,
        output_data,
        configurations,
        dΩ;
        tasks_per_thread = tasks_per_thread,
        observer_four_velocity = observer_four_velocity,
        flux_direction = flux_direction)
    return nothing
end

function fluxes_task(chunk,
    F::AbstractVector,
    initial_data,
    output_data,
    configurations,
    dΩ;
    kwargs...)
    cache = fluxes_task_init(configurations; kwargs...)
    for i in chunk
        @views begin
            ki = initial_data[5:8, i]
        end
        nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
        nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
        F[i] *= nu * nn * dΩ[i]
    end
    return nothing
end

function fluxes_task(chunk,
    F::AbstractMatrix,
    initial_data,
    output_data,
    configurations,
    dΩ;
    kwargs...)
    cache = fluxes_task_init(configurations; kwargs...)
    for i in chunk
        @views begin
            ki = initial_data[5:8, i]
        end
        nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
        nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
        F[:, i] *= nu * nn * dΩ[i]
    end
    return nothing
end

function fluxes_checks(F::AbstractArray,
    initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations)
    nrays = size(F, ndims(F))
    eight_components(initial_data, output_data) ||
        throw(DimensionMismatch("The initial and output data must have eight components."))
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))
    size(initial_data, 2) == nrays ||
        throw(DimensionMismatch("The number of rays in the initial and output data must be the same as the length of F."))
    number_of_pixels(configurations.camera) == nrays ||
        throw(DimensionMismatch("The number of pixels in the camera must be the same as the length of F."))
    return nothing
end

function fluxes_task_init(configurations; observer_four_velocity, flux_direction)
    spacetime = configurations.spacetime
    cache = postprocess_cache(configurations)
    camera = configurations.camera
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity)
    flux_direction!(cache, flux_direction, camera, spacetime)
    return cache
end
