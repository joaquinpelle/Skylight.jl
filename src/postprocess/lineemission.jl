"""
    line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; 
                           emission_profile::Function, bin_size::Number=NaN, num_bins::Int=NaN,
                           start::Number=NaN, stop::Number=NaN)

Compute the binned intensity of a line emission spectrum.

# Arguments
- `initial_data`: Initial condition data.
- `output_data`: Output data from the radiation model.
- `configurations::VacuumOTEConfigurations`: Configuration parameters for the model.

# Keywords
- `emission_profile::Function`: User-defined function describing the emission profile.
- `bin_size::Union{Nothing, Number}=nothing`: Size of each bin. Either `bin_size` or `num_bins` must be specified.
- `num_bins::Union{Nothing, Number}=nothing`: Number of bins. Either `bin_size` or `num_bins` must be specified.
- `start::Union{Nothing, Number}=nothing`: Lower bound of the range to be binned. If unspecified, the minimum of the energy quotients will be used.
- `stop::Union{Nothing, Number}=nothing`: Upper bound of the range to be binned. If unspecified, the maximum of the energy quotients will be used.
- `observer_four_velocity::AbstractVector` (optional): The four-velocity of the observer. If not provided, a default static four-velocity will be used.
- `flux_direction::AbstractVector` (optional): The direction in which to measure the flu. If not provided, a default direction will be used.

# Returns
- `binned_fluxes`: Array of the binned intensity in each bin.
- `bins`: Array of the bin edges.

# Notes
`observer_four_velocity` and `flux_direction` are only accepted if `configurations.camera` is a `PinholeCamera`.
"""
line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; kwargs...) = line_emission_spectrum(initial_data, output_data, configurations, configurations.camera; kwargs...)

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

    not_simultaneously_nothing(num_bins, bin_size) || throw(ArgumentError("Either bin_size or num_bins must be specified."))
    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(configurations)

    Nrays = size(initial_data, 2)
    F = zeros(Nrays)
    q = zeros(Nrays)
    at_source = zeros(Bool, Nrays)

    for i in axes(initial_data, 2)

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

    not_simultaneously_nothing(num_bins, bin_size) || throw(ArgumentError("Either bin_size or num_bins must be specified."))
    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(configurations)

    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    flux_direction!(cache, flux_direction, camera, spacetime) 

    dΩ = pixel_solid_angles(camera)
    Nrays = size(initial_data, 2)
    F = zeros(Nrays)
    q = zeros(Nrays)
    at_source = zeros(Bool, Nrays)

    for i in axes(initial_data, 2)

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