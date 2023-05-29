"""
    line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; 
                           emission_profile::Function, bin_size::Number=NaN, num_bins::Int=NaN,
                           start::Number=NaN, stop::Number=NaN)

Compute the binned intensity of a line emission spectrum.

# Arguments
- `initial_data`: Initial condition data.
- `output_data`: Output data from the radiation model.
- `configurations::VacuumOTEConfigurations`: Configuration parameters for the model.
- `emission_profile::Function`: User-defined function describing the emission profile.

# Keywords
- `bin_size::Union{Nothing, Number}=nothing`: Size of each bin. Either `bin_size` or `num_bins` must be specified.
- `num_bins::Union{Nothing, Number}=nothing`: Number of bins. Either `bin_size` or `num_bins` must be specified.
- `start::Union{Nothing, Number}=nothing`: Lower bound of the range to be binned. If unspecified, the minimum of the energy quotients will be used.
- `stop::Union{Nothing, Number}=nothing`: Upper bound of the range to be binned. If unspecified, the maximum of the energy quotients will be used.
# Returns
- `binned_fluxes`: Array of the binned intensity in each bin.
- `bins`: Array of the bin edges.
"""
function line_emission_spectrum(
    initial_data, 
    output_data, 
    configurations::VacuumOTEConfigurations; 
    emission_profile::Function, 
    bin_size::Union{Number,Nothing}=nothing, 
    num_bins::Union{Int,Nothing}=nothing,
    start::Union{Number,Nothing}=nothing, 
    stop::Union{Number,Nothing}=nothing)

    if bin_size===nothing && num_bins===nothing 
        throw(ArgumentError("Either bin_size or num_bins must be specified."))
    end

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

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

        dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
        F[i] = q[i]^3*emission_profile(pf, spacetime, model)
    end

    #If start or stop not provided, set to extremes
    if start===nothing
        start = minimum(q[at_source])
    end
    if stop===nothing
        stop = maximum(q[at_source])
    end
    
    num_bins = Int(round((stop-start)/bin_size))

    bins = create_bins(bin_size=bin_size, num_bins=num_bins, start=start, stop=stop)
    binned_fluxes = bin_values_and_sum_weights(bins, q[at_source], F[at_source])
    
    normalize_by_pixel_area!(binned_fluxes, configurations)
    normalize_by_image_plane_distance!(binned_fluxes, configurations)
    return binned_fluxes, bins
end

"""
    line_emission_spectrum(initial_data, output_data, configurations::VacuumOTEConfigurations; 
                           emission_profile::Function, bin_size::Number=NaN, num_bins::Int=NaN,
                           start::Number=NaN, stop::Number=NaN)

Compute the binned intensity of a line emission spectrum.

# Arguments
- `initial_data`: Initial condition data.
- `output_data`: Output data from the radiation model.
- `configurations::VacuumOTEConfigurations`: Configuration parameters for the model.
- `emission_profile::Function`: User-defined function describing the emission profile.

# Keywords
- `start::Union{Number, Nothing}`: Lower bound of the range to be binned. If unspecified, the minimum of the energy quotients will be used.
- `stop::Union{Nothing, Number}=nothing`: Upper bound of the range to be binned. If unspecified, the maximum of the energy quotients will be used.
- `bin_size_conditioner::Number`: Conditioner of the bin size. The bin size is required to be larger than the conditioner times the max variation of the energy factor.
- `edge_width::Int=3`: Width of the edge to be ignored for bin size conditioning (since edges have unusually large local variations, especially in the presence of light-rings and similar phenomena).
# Returns
- `binned_fluxes`: Array of the binned intensity in each bin.
- `bins`: Array of the bin edges.
"""
function line_emission_spectrum(
    initial_data, 
    output_data, 
    configurations::VacuumOTEConfigurations; 
    emission_profile::Function, 
    start::Union{Number,Nothing}=nothing, 
    stop::Union{Number,Nothing}=nothing,
    bin_size_conditioner::Number,
    edge_width::Int)

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

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

        dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
        F[i] = q[i]^3*emission_profile(pf, spacetime, model)
    end

    #If start or stop not provided, set to extremes
    if start===nothing
        start = minimum(q[at_source])
    end
    if stop===nothing
        stop = maximum(q[at_source])
    end

    num_bins = infer_num_bins(q, at_source, start, stop, bin_size_conditioner, edge_width, configurations.image_plane)

    bins = create_bins(bin_size=bin_size, num_bins=num_bins, start=start, stop=stop)
    binned_fluxes = bin_values_and_sum_weights(bins, q[at_source], F[at_source])
    
    normalize_by_pixel_area!(binned_fluxes, configurations)
    normalize_by_image_plane_distance!(binned_fluxes, configurations)
    return binned_fluxes, bins
end

"""
    infer_num_bins(q, at_source, start, stop, bin_size_conditioner, edge_width, image_plane)

Infer the number of bins from the energy quotients and the bin size conditioner. It chooses the larges number of bins that satisfies
that the bin size is larger than the conditioner times the maximum local variation of the energy quotients, for a given approximation of the local variation.

# Arguments

- `q`: Array of energy quotients.
- `at_source`: Array of booleans indicating whether the final position is at the source.
- `start`: Lower bound of the range to be binned.
- `stop`: Upper bound of the range to be binned.
- `bin_size_conditioner`: Conditioner of the bin size.
- `edge_width`: Width of the edge to be ignored for bin size conditioning (since edges have unusually large local variations, especially in the presence of light-rings and similar phenomena).
- `image_plane`: Image plane.

# Returns

- `num_bins`: Number of bins.
"""
function infer_num_bins(q, at_source, start, stop, bin_size_conditioner, edge_width, image_plane)
    
    Nα, Nβ = numbers_of_nodes_per_side(image_plane)
    dα, dβ = grid_spacing(image_plane)
    
    at_edge = detect_edges(edge_width, reshape(at_source, Nα, Nβ))
    
    δq = approximate_gradient_norm(reshape(q, Nα, Nβ), dα, dβ)
    max_δq = maximum(δq[.!at_edge])
    return Int(floor((stop-start)/(bin_size_conditioner*max_δq))) 
end

