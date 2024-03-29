"""
    line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; 
                           bin_size::Number=NaN, num_bins::Int=NaN,
                           start::Number=NaN, stop::Number=NaN)

Compute the binned intensity of a line emission spectrum.

# Arguments
- `initial_data`: Initial condition data.
- `output_data`: Output data from the radiation model.
- `configurations::VacuumOTEConfigurations`: Configuration parameters for the model.

# Keywords
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
function line_emission_spectrum(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations;
    kwargs...)
    line_emission_spectrum(initial_data,
        output_data,
        configurations,
        configurations.camera;
        kwargs...)
end

function line_emission_spectrum(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    ::ImagePlane;
    num_bins::Union{Number, Nothing} = nothing,
    bin_size::Union{Number, Nothing} = nothing,
    start::Union{Number, Nothing} = nothing,
    stop::Union{Number, Nothing} = nothing,
    tasks_per_thread::Int = 2)
    not_simultaneously_nothing(num_bins, bin_size) ||
        throw(ArgumentError("Either bin_size or num_bins must be specified."))
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    nrays = size(initial_data, 2)
    F = zeros(nrays)
    q = zeros(nrays)

    at_source = zeros(Bool, nrays)
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
                if !is_final_position_at_source(pf, spacetime, model)
                    continue
                end
                at_source[i] = true
                metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
                q[i] = energies_quotient(ki, kf, cache) #TODO add energy as argument below
                F[i] = q[i]^3 * line_emission_profile(pf,
                    -kf,
                    cache.rest_frame_four_velocity,
                    cache.emitter_metric,
                    spacetime,
                    model,
                    coords_top,
                    cache)
            end
        end
    end
    if start === nothing
        start = minimum(q[at_source])
    end
    if stop === nothing
        stop = maximum(q[at_source])
    end
    bins_edges = create_bins(bin_size = bin_size, num_bins = num_bins, start = start, stop = stop)
    binned_fluxes = bin_values_and_sum_weights(edges=bins_edges, values=q[at_source], weights=F[at_source])
    normalize_by_pixel_area!(binned_fluxes, configurations)
    normalize_by_camera_distance!(binned_fluxes, configurations)
    return binned_fluxes, bins_edges
end

function line_emission_spectrum(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    camera::PinholeCamera;
    num_bins::Union{Number, Nothing} = nothing,
    bin_size::Union{Number, Nothing} = nothing,
    start::Union{Number, Nothing} = nothing,
    stop::Union{Number, Nothing} = nothing,
    observer_four_velocity::Union{AbstractVector, Nothing} = nothing,
    flux_direction::Union{AbstractVector, Nothing} = nothing,
    tasks_per_thread::Int = 2)
    not_simultaneously_nothing(num_bins, bin_size) ||
        throw(ArgumentError("Either bin_size or num_bins must be specified."))
    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("initial_data and output_data must have the same size."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    dΩ = pixel_solid_angles(camera)
    nrays = size(initial_data, 2)
    F = zeros(nrays)
    q = zeros(nrays)

    at_source = zeros(Bool, nrays)

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
                at_source[i] = true
                emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
                nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
                nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
                q[i] = energies_quotient(ki, kf, cache)
                F[i] = nu * nn * q[i]^3 * dΩ[i] *
                       line_emission_profile(pf,
                           -kf,
                           cache.rest_frame_four_velocity,
                           cache.emitter_metric,
                           spacetime,
                           model,
                           coords_top,
                           cache)
            end
        end
    end
    if start === nothing
        start = minimum(q[at_source])
    end
    if stop === nothing
        stop = maximum(q[at_source])
    end
    bins_edges = create_bins(bin_size = bin_size, num_bins = num_bins, start = start, stop = stop)
    binned_fluxes = bin_values_and_sum_weights(edges=bins_edges, values=q[at_source], weights=F[at_source])
    return binned_fluxes, bins_edges
end
