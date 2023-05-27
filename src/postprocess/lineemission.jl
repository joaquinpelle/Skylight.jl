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
- `bin_size::Number=NaN`: Size of each bin. Either `bin_size` or `num_bins` must be specified.
- `num_bins::Int=NaN`: Number of bins. Either `bin_size` or `num_bins` must be specified.
- `start::Number=NaN`: Lower bound of the range to be binned. If unspecified, the minimum of the energy quotients will be used.
- `stop::Number=NaN`: Upper bound of the range to be binned. If unspecified, the maximum of the energy quotients will be used.

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

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    observed_specific_fluxes = zeros(Nrays)
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
        
        q[i] = get_energies_quotient(ki, kf, cache)

        observed_specific_fluxes[i] = q[i]^3*emission_profile(pf, spacetime, model)

    end

    q = q[at_source]
    observed_specific_fluxes = observed_specific_fluxes[at_source] 
    
    #If start or stop not provided, set to extremes
    start===nothing && (start = minimum(q)) 
    stop===nothing && (stop = maximum(q))  
    
    bins = create_bins(bin_size=bin_size, num_bins=num_bins, start=start, stop=stop)

    binned_fluxes = bin_values_and_sum_weights(bins, q, observed_specific_fluxes)
    
    normalize_by_pixel_area!(binned_fluxes, configurations)
    normalize_by_image_plane_distance!(binned_fluxes, configurations)
    
    return binned_fluxes, bins
end