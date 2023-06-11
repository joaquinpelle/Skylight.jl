observed_bolometric_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations; kwargs...) = observed_bolometric_intensities(initial_data, output_data, configurations, configurations.camera; kwargs...)
observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations, energy::Number; kwargs...) = observed_specific_intensities(initial_data, output_data, configurations, [energy]; kwargs...)
observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations, energies::AbstractVector; kwargs...) = observed_specific_intensities(initial_data, output_data, configurations, configurations.camera, energies; kwargs...)

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
function observed_bolometric_intensities(initial_data::AbstractMatrix, output_data::AbstractMatrix, configurations::VacuumOTEConfigurations, ::ImagePlane)

    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have eight components.")) 
    
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(configurations)

    Nrays = size(initial_data, 2)
    q = zeros(Nrays)
    Iobs = zeros(Nrays)

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
        
        metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)

        #The difference with the ETO scheme here should be the minus sign in front of the final momentum
        #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
        Iobs[i] = q[i]^4*emitted_bolometric_intensity(pf, -kf, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
    end
    return Iobs, q
end

""" observed_bolometric_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations, camera::PinholeCamera; observer_four_velocity=nothing) 

Compute observed bolometric intensities and energy quotients for a set of rays defined by the initial and final conditions. The function also checks whether the final position of each ray is at the source, in which case it proceeds with computation, else it skips to the next ray and the values are set to zero. 

    # Arguments
- `initial_data::AbstractMatrix`: A matrix containing the initial data of the rays. The first four entries for each ray represent the initial position, while entries five to eight represent the initial momentum.
- `output_data::AbstractMatrix`: A matrix containing the final data of the rays. The first four entries for each ray represent the final position, while entries five to eight represent the final momentum.
- `configurations::VacuumOTEConfigurations`: The configurations with which the initial and output data were obtained.
- `camera::PinholeCamera`: The camera model used to capture the scene.
- `observer_four_velocity::AbstractVector` (optional): The four-velocity of the observer. If not provided, a default static four-velocity will be used.

# Returns
- `Iobs::Vector`: A vector of observed bolometric intensities for each ray, normalized by the distance to the image plane. The intensities are in CGS units.
- `q::Vector`: A vector of energy quotients for each ray, representing the ratio of final energy in the local emission frame to initial energy in the observer frame.

# Notes
    Output units are CGS. The observer four-velocity and flux direction, if provided, must satisfy the conditions of being timelike and spacelike, respectively, as per the observer metric.
"""
function observed_bolometric_intensities(initial_data::AbstractMatrix, output_data::AbstractMatrix, configurations::VacuumOTEConfigurations, camera::PinholeCamera; observer_four_velocity=nothing)
    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have eight components.")) 

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 

    Nrays = size(initial_data, 2) 
    q = zeros(Nrays)
    Iobs = zeros(Nrays)

    for i in axes(initial_data, 2)

        @views begin 
            ki = initial_data[5:8,i]
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]
        end

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end

        emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
        
        #The difference with the ETO scheme here should be the minus sign in front of the final momentum
        #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
        Iem = emitted_bolometric_intensity(pf, -kf, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
        Iobs[i] = q[i]^4*Iem        
    end
    return Iobs, q
end

"""
    observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations)

Compute observed specific intensities and energy quotients for a set of rays defined by the initial and final conditions. The function also checks whether the final position of each ray is at the source, in which case it proceeds with computation, else it skips to the next ray and the values are set to zero.

# Arguments
- `initial_data::AbstractMatrix`: A matrix containing the initial data of the rays. The first four entries for each ray represent the initial position, while entries five to eight represent the initial momentum.
- `output_data::AbstractMatrix`: A matrix containing the final data of the rays. The first four entries for each ray represent the final position, while entries five to eight represent the final momentum.
- `configurations::VacuumOTEConfigurations`: The configurations with which the initial and output data were obtained.
- `observation_energies::Vector`: A vector of energies in CGS units at which the specific intensities are to be computed.

# Returns
- `Iobs::Vector`: A vector of observed specific intensities in CGS units for each ray, normalized by the distance to the image plane. The intensities are in CGS units.
- `q::Vector`: A vector of energy quotients for each ray, representing the ratio of final energy in the local emisison frame to initial energy in the observer frame.

# Notes
    Input energy units must be CGS. Output units are CGS.
"""
function observed_specific_intensities(initial_data::AbstractMatrix, output_data::AbstractMatrix, configurations::VacuumOTEConfigurations, ::ImagePlane, observation_energies::AbstractVector)
    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have eight components.")) 
    
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(configurations)

    Nrays = size(initial_data, 2)
    NE = length(observation_energies)

    q = zeros(Nrays)
    Iobs = zeros(NE, Nrays)
    for i in axes(initial_data,2)

        @views begin 
            pi = initial_data[1:4,i]
            ki = initial_data[5:8,i]
            
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]
        end

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end

        metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
        
        for j in axes(observation_energies, 1)
            emitted_energy = observation_energies[j]/q[i]

            #The difference with the ETO scheme here should be the minus sign in front of the final momentum
            #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
            Iobs[j, i] = q[i]^3*emitted_specific_intensity(pf, -kf, emitted_energy, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
        end
    end
    return Iobs, q
end

"""
    observed_specific_intensities(initial_data, output_data, configurations::VacuumOTEConfigurations, observation_energies; observer_four_velocity=nothing)

Compute observed specific intensities and energy quotients for a set of rays defined by the initial and final conditions. The function also checks whether the final position of each ray is at the source, in which case it proceeds with computation, else it skips to the next ray and the values are set to zero.

# Arguments
- `initial_data::AbstractMatrix`: A matrix containing the initial data of the rays. The first four entries for each ray represent the initial position, while entries five to eight represent the initial momentum.
- `output_data::AbstractMatrix`: A matrix containing the final data of the rays. The first four entries for each ray represent the final position, while entries five to eight represent the final momentum.
- `configurations::VacuumOTEConfigurations`: The configurations with which the initial and output data were obtained.
- `observation_energies::Vector`: A vector of energies in CGS units at which the specific intensities are to be computed.
- `observer_four_velocity::AbstractVector` (optional): The four-velocity of the observer. If not provided, a default static four-velocity will be used.

# Returns
- `Iobs::Matrix`: A matrix of observed specific intensities in CGS units for each ray, normalized by the distance to the image plane. The intensities are in CGS units.
- `q::Vector`: A vector of energy quotients for each ray, representing the ratio of final energy in the local emission frame to initial energy in the observer frame.

# Notes
    Input energy units must be CGS. Output units are CGS. The observer four-velocity and flux direction, if provided, must satisfy the conditions of being timelike and spacelike, respectively, as per the observer metric.
"""
function observed_specific_intensities(initial_data::AbstractMatrix, output_data::AbstractMatrix, configurations::VacuumOTEConfigurations, camera::PinholeCamera, observation_energies; observer_four_velocity=nothing)
    
    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have eight components."))
    
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)

    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 

    Nrays = size(initial_data, 2)
    NE = length(observation_energies)

    q = zeros(Nrays)
    Iobs = zeros(NE, Nrays)
    for i in axes(initial_data, 2)

        @views begin 
            ki = initial_data[5:8,i]
            pf = output_data[1:4,i]
            kf = output_data[5:8,i]
        end

        if !is_final_position_at_source(pf, spacetime, model)
            continue
        end

        emitter_metric_and_four_velocity!(cache, pf, spacetime, model, coords_top)
        q[i] = energies_quotient(ki, kf, cache)
        
        for j in axes(observation_energies, 1)
            emitted_energy = observation_energies[j]/q[i]
            #The difference with the ETO scheme here should be the minus sign in front of the final momentum
            #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
            Iem = emitted_specific_intensity(pf, -kf, emitted_energy, cache.emitter_four_velocity, cache.emitter_metric, spacetime, model, coords_top)
            Iobs[j, i] = q[i]^3*Iem            
        end
    end
    return Iobs, q
end