function line_emission_spectrum(
    initial_data::AbstractMatrix, 
    output_data::AbstractMatrix, 
    configurations::VacuumOTEConfigurations,
    ::ImagePlane; 
    emission_profile::Function, 
    num_bins::Union{Number,Nothing}=nothing,
    bin_size::Union{Number,Nothing}=nothing, 
    start::Union{Number,Nothing}=nothing, 
    stop::Union{Number,Nothing}=nothing)

    if bin_size===nothing && num_bins===nothing 
        throw(ArgumentError("Either bin_size or num_bins must be specified."))
    end

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(configurations)

    Nrays = number_of_initial_conditions(configurations)
    F = zeros(Nrays)
    q = zeros(Nrays)
    at_source = zeros(Bool, Nrays)

    for i in 1:Nrays

        @views begin 

            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]

        end

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end
        at_source[i] = true

        metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
        F[i] = q[i]^3*emission_profile(pf, spacetime, model)
    end

    if start===nothing
        start = minimum(q[at_source])
    end
    if stop===nothing
        stop = maximum(q[at_source])
    end
    
    bins = create_bins(bin_size=bin_size, num_bins=num_bins, start=start, stop=stop)
    binned_fluxes = bin_values_and_sum_weights(bins, q[at_source], F[at_source])
    
    normalize_by_pixel_area!(binned_fluxes, configurations)
    normalize_by_image_plane_distance!(binned_fluxes, configurations)
    return binned_fluxes, bins
end

function line_emission_spectrum(
    initial_data::AbstractMatrix, 
    output_data::AbstractMatrix, 
    configurations::VacuumOTEConfigurations,
    camera::PinholeCamera; 
    emission_profile::Function, 
    num_bins::Union{Number,Nothing}=nothing,
    bin_size::Union{Number,Nothing}=nothing, 
    start::Union{Number,Nothing}=nothing, 
    stop::Union{Number,Nothing}=nothing,
    observer_four_velocity::AbstractVector=nothing,
    flux_direction::AbstractVector=nothing)

    if bin_size===nothing && num_bins===nothing 
        throw(ArgumentError("Either bin_size or num_bins must be specified."))
    end

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(configurations)

    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    flux_direction!(cache, flux_direction, camera, spacetime) 

    dΩ = pixel_solid_angles(camera)
    Nrays = number_of_initial_conditions(configurations)
    F = zeros(Nrays)
    q = zeros(Nrays)
    at_source = zeros(Bool, Nrays)

    for i in 1:Nrays

        @views begin 
            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]
        end

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end
        at_source[i] = true

        emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
        nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
        nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
        q[i] = energies_quotient(ki, kf, cache)
        F[i] = nu*nn*q[i]^3*emission_profile(pf, spacetime, model)*dΩ[i]
    end

    if start===nothing
        start = minimum(q[at_source])
    end
    if stop===nothing
        stop = maximum(q[at_source])
    end
    
    bins = create_bins(bin_size=bin_size, num_bins=num_bins, start=start, stop=stop)
    binned_fluxes = bin_values_and_sum_weights(bins, q[at_source], F[at_source])
    return binned_fluxes, bins
end