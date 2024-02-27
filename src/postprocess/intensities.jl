function observed_bolometric_intensities(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations;
    kwargs...)
    observed_bolometric_intensities(initial_data,
        output_data,
        configurations,
        configurations.camera;
        kwargs...)
end
function observed_specific_intensities(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    energy::Real;
    kwargs...)
    observed_specific_intensities(initial_data,
        output_data,
        configurations,
        [energy];
        kwargs...)
end
function observed_specific_intensities(initial_data,
    output_data,
    configurations::VacuumOTEConfigurations,
    energies::AbstractVector;
    kwargs...)
    observed_specific_intensities(initial_data,
        output_data,
        configurations,
        configurations.camera,
        energies;
        kwargs...)
end

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

function observed_bolometric_intensities(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    ::ImagePlane;
    tasks_per_thread::Int = 2)
    Iobs = postprocess_init(initial_data, output_data, configurations)
    function task(chunk, Iobs, initial_data, output_data, configurations)
        spacetime, model, coords_top, cache = task_init(configurations)
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
            metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
            q = energies_quotient(ki, kf, cache)
            #The difference with the ETO scheme here should be the minus sign in front of the final momentum
            #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
            Iobs[i] = q^4 * rest_frame_bolometric_intensity(pf,
                -kf,
                cache.rest_frame_four_velocity,
                cache.emitter_metric,
                spacetime,
                model,
                coords_top,
                cache.model_cache)
        end
        return nothing
    end
    itr = axes(Iobs, ndims(Iobs))
    tmap(task,
        itr,
        Iobs,
        initial_data,
        output_data,
        configurations;
        tasks_per_thread = tasks_per_thread)
    return Iobs
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
function observed_bolometric_intensities(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    ::PinholeCamera;
    observer_four_velocity = nothing,
    tasks_per_thread::Int = 2)
    Iobs = postprocess_init(initial_data, output_data, configurations)
    function task(chunk,
        Iobs,
        initial_data,
        output_data,
        configurations;
        observer_four_velocity)
        spacetime, model, coords_top, cache = task_init(configurations;
            observer_four_velocity)
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
            #The difference with the ETO scheme here should be the minus sign in front of the final momentum
            #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
            Iem = rest_frame_bolometric_intensity(pf,
                -kf,
                cache.rest_frame_four_velocity,
                cache.emitter_metric,
                spacetime,
                model,
                coords_top,
                cache.model_cache)
            Iobs[i] = q^4 * Iem
        end
        return nothing
    end
    itr = axes(Iobs, ndims(Iobs))
    tmap(task,
        itr,
        Iobs,
        initial_data,
        output_data,
        configurations;
        observer_four_velocity = observer_four_velocity,
        tasks_per_thread = tasks_per_thread)
    return Iobs
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
function observed_specific_intensities(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    ::ImagePlane,
    observation_energies::AbstractVector;
    tasks_per_thread::Int = 2)
    Iobs = postprocess_init(initial_data, output_data, configurations, observation_energies)
    function task(chunk,
        Iobs,
        initial_data,
        output_data,
        configurations,
        observation_energies)
        spacetime, model, coords_top, cache = task_init(configurations)
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
            metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
            q = energies_quotient(ki, kf, cache)
            for j in axes(observation_energies, 1)
                emitted_energy = observation_energies[j] / q
                #The difference with the ETO scheme here should be the minus sign in front of the final momentum
                #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
                Iobs[j, i] = q^3 * rest_frame_specific_intensity(pf,
                    -kf,
                    emitted_energy,
                    cache.rest_frame_four_velocity,
                    cache.emitter_metric,
                    spacetime,
                    model,
                    coords_top,
                    cache.model_cache)
            end
        end
    end
    itr = axes(Iobs, ndims(Iobs))
    tmap(task,
        itr,
        Iobs,
        initial_data,
        output_data,
        configurations,
        observation_energies;
        tasks_per_thread = tasks_per_thread)
    return Iobs
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
function observed_specific_intensities(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    ::PinholeCamera,
    observation_energies;
    observer_four_velocity = nothing,
    tasks_per_thread::Int = 2)
    Iobs = postprocess_init(initial_data, output_data, configurations, observation_energies)
    function task(chunk,
        Iobs,
        initial_data,
        output_data,
        configurations,
        observation_energies;
        observer_four_velocity)
        spacetime, model, coords_top, cache = task_init(configurations;
            observer_four_velocity)
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
            for j in axes(observation_energies, 1)
                emitted_energy = observation_energies[j] / q
                #The difference with the ETO scheme here should be the minus sign in front of the final momentum
                #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
                Iem = rest_frame_specific_intensity(pf,
                    -kf,
                    emitted_energy,
                    cache.rest_frame_four_velocity,
                    cache.emitter_metric,
                    spacetime,
                    model,
                    coords_top,
                    cache.model_cache)
                Iobs[j, i] = q^3 * Iem
            end
        end
    end
    itr = axes(Iobs, ndims(Iobs))
    tmap(task,
        itr,
        Iobs,
        initial_data,
        output_data,
        configurations,
        observation_energies;
        observer_four_velocity = observer_four_velocity,
        tasks_per_thread = tasks_per_thread)
    return Iobs
end

"""
    observed_specific_intensities(initial_data, output_data, configurations::NonVacuumOTEConfigurations)

Compute observed specific intensities

# Arguments
- `output_data::AbstractMatrix`: A matrix containing the final data of the rays. The first four entries for each ray represent the final position, while entries five to eight represent the final momentum.
- `configurations::NonVacuumOTEConfigurations`: The configurations with which the initial and output data were obtained.

# Returns
- `Iobs::AbstractMatrix`: A vector of observed specific intensities in CGS units for each ray, normalized by the distance to the image plane. The intensities are in CGS units.

"""
#TODO add itr axes(Iobs)
@kwdispatch observed_specific_intensities(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::NonVacuumOTEConfigurations)

@kwmethod function observed_specific_intensities(::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::NonVacuumOTEConfigurations;)
    NE = length(configurations.observation_energies)
    return @. observation_energies^3 * output_data[(9 + NE):end, :]
end

@kwmethod function observed_specific_intensities(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::NonVacuumOTEConfigurations;
    observer_four_velocity)
    Iobs = postprocess_init(initial_data, output_data, configurations)
    function task(chunk,
        Iobs,
        initial_data,
        output_data,
        configurations;
        observer_four_velocity)
        observer_metric = metric(configurations.camera.position, configurations.spacetime)
        NE = length(configurations.observation_energies)
        for i in chunk
            @views begin
                ki = initial_data[5:8, i]
            end
            observer_rest_frame_energy = scalar_product(ki,
                observer_four_velocity,
                observer_metric)
            @. Iobs[:, i] = (observation_energies * observer_rest_frame_energy)^3 *
                            output_data[(9 + NE):end, i]
        end
    end
    tmap(task,
        axes(Iobs, ndims(Iobs)),
        Iobs,
        initial_data,
        output_data,
        configurations;
        observer_four_velocity = observer_four_velocity,
        tasks_per_thread = tasks_per_thread)
    return Iobs
end

function observed_bolometric_intensities_serial(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumOTEConfigurations,
    ::ImagePlane)
    Iobs = postprocess_init(initial_data, output_data, configurations)
    spacetime, model, coords_top, cache = task_init(configurations)
    @inbounds begin
        for i in axes(initial_data, 2)
            @views begin
                pi = initial_data[1:4, i]
                ki = initial_data[5:8, i]
                pf = output_data[1:4, i]
                kf = output_data[5:8, i]
            end
            if !is_final_position_at_source(pf, spacetime, model)
                continue
            end
            metrics_and_four_velocities!(cache, pi, pf, spacetime, model, coords_top)
            q = energies_quotient(ki, kf, cache)
            #The difference with the ETO scheme here should be the minus sign in front of the final momentum
            #at get emitted intensity, and the is_final_position_at_source call (at observer in ETO)...
            Iobs[i] = q^4 * rest_frame_bolometric_intensity(pf,
                -kf,
                cache.rest_frame_four_velocity,
                cache.emitter_metric,
                spacetime,
                model,
                coords_top,
                cache.model_cache)
        end
    end
    return Iobs
end
