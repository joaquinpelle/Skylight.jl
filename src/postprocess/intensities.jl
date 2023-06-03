"""
    observed_bolometric_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations)

Compute observed bolometric intensities and energy quotients for a set of rays defined by the initial and final conditions. The function also checks whether the final position of each ray is at the source, in which case it proceeds with computation, else it skips to the next ray and the values are set to zero.

# Arguments
- `initial_data::AbstractMatrix`: A matrix containing the initial data of the rays. The first four entries for each ray represent the initial position, while entries five to eight represent the initial momentum.
- `output_data::AbstractMatrix`: A matrix containing the final data of the rays. The first four entries for each ray represent the final position, while entries five to eight represent the final momentum.
- `configurations::VacuumOTEConfigurations`: The configurations with which the initial and output data were obtained.

# Returns
- `Iobs::Vector`: A vector of observed bolometric intensities for each ray, normalized by the distance to the image plane. The intensities are in CGS units.
- `q::Vector`: A vector of energy quotients for each ray, representing the ratio of final energy in the local emisison frame to initial energy in the observer frame.

# Notes
    Output units are CGS.
"""
function observed_bolometric_intensities(initial_data::AbstractMatrix, output_data::AbstractMatrix, configurations::VacuumOTEConfigurations)

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    q = zeros(Nrays)
    Iobs = zeros(Nrays)

    @inbounds begin
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

            dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coords_top)
            q[i] = energies_quotient(ki, kf, cache)

            #The difference with the ETO scheme here should be the minus sign in front of the final momentum
            #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
            Iobs[i] = q[i]^4*emitted_bolometric_intensity(pf, -kf, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
        end
    end

    normalize_by_image_plane_distance!(Iobs, configurations)
    return Iobs, q
end

"""
    observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations)

Compute observed specific intensities and energy quotients for a set of rays defined by the initial and final conditions. The function also checks whether the final position of each ray is at the source, in which case it proceeds with computation, else it skips to the next ray and the values are set to zero.

# Arguments
- `initial_data::AbstractMatrix`: A matrix containing the initial data of the rays. The first four entries for each ray represent the initial position, while entries five to eight represent the initial momentum.
- `output_data::AbstractMatrix`: A matrix containing the final data of the rays. The first four entries for each ray represent the final position, while entries five to eight represent the final momentum.
- `configurations::VacuumOTEConfigurations`: The configurations with which the initial and output data were obtained.
- `observed_energies::Vector`: A vector of energies in CGS units at which the specific intensities are to be computed.

# Returns
- `Iobs::Vector`: A vector of observed specific intensities in CGS units for each ray, normalized by the distance to the image plane. The intensities are in CGS units.
- `q::Vector`: A vector of energy quotients for each ray, representing the ratio of final energy in the local emisison frame to initial energy in the observer frame.

# Notes
    Input energy units must be CGS. Output units are CGS.
"""
function observed_specific_intensities(initial_data::AbstractMatrix, output_data::AbstractMatrix, configurations::VacuumOTEConfigurations, observed_energies)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = get_postprocess_cache(configurations)
    dump_observer_four_velocity_in!(cache)

    Nrays = number_of_initial_conditions(configurations)
    NE = length(observed_energies)

    q = zeros(Nrays)
    Iobs = zeros(NE, Nrays)
    @inbounds begin
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

            dump_metrics_and_emitter_four_velocity_in!(cache, pi, pf, spacetime, model, coords_top)
            q[i] = energies_quotient(ki, kf, cache)
            
            for j in 1:NE
                emitted_energy = observed_energies[j]/q[i]

                #The difference with the ETO scheme here should be the minus sign in front of the final momentum
                #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
                Iobs[j, i] = q[i]^3*emitted_specific_intensity(pf, -kf, emitted_energy, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
            end
        end
    end

    normalize_by_image_plane_distance!(Iobs, configurations)
    return Iobs, q
end