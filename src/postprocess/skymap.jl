function specific_flux_skymap(initial_data::AbstractMatrix,
    output_data::AbstractMatrix,
    configurations::VacuumETOConfigurations,
    observation_energies::AbstractVector; #must be in CGS
    Nθ::Int,
    Nϕ::Int,
    tasks_per_thread::Int = 2)

    same_size(initial_data, output_data) ||
        throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) ||
        throw(DimensionMismatch("The initial and output data must have eight components."))

    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    period = system_period(model)
    θbins = create_bins(num_bins = Nθ, start = 0.0, stop = π)
    ϕbins = create_bins(num_bins = Nϕ, start = 0.0, stop = 1.0)

    itr = axes(initial_data, ndims(initial_data))
    # Break the work into chunks. More chunks per thread has better load balancing but more overhead
    chunk_size = max(1, length(itr) ÷ (tasks_per_thread * nthreads()))
    # Map over the chunks, creating an array of spawned tasks. Sync to wait for the tasks to finish.
    tasks = map(Iterators.partition(itr, chunk_size)) do chunk
        @spawn begin
            # Initialize an array to hold the sum of `q` values in each bin
            Fsums = zeros(length(observation_energies), length(θbins)-1, length(ϕbins)-1)
            cache = postprocess_cache(configurations)
            for i in chunk
                @views begin
                    pi = initial_data[1:4, i]
                    ki = initial_data[5:8, i]
                    pf = output_data[1:4, i]
                    kf = output_data[5:8, i]
                end
                if !is_final_position_at_observer(pf, configurations)
                    continue
                end
                emitter_metric_four_velocity_and_normal!(cache, pi, spacetime, model, coords_top)
                #Assuming the observation frame is static, and the photon packages are emitted with unit frequency in the rest frame
                q = kf[1] 
                ϕ = observed_phase(pf, period, coords_top)
                θbin_index = searchsortedlast(θbins, pf[3])
                ϕbin_index = searchsortedlast(ϕbins, ϕ)
                for j in axes(observation_energies, 1)
                    Eobs = observation_energies[j]
                    Eem = Eobs/q
                    #TODO check units
                    Fvalue = Eobs*photon_package_weight(pi, 
                        ki, 
                        Eem, 
                        cache.emitter_metric, 
                        cache.rest_frame_four_velocity,
                        cache.surface_normal,
                        spacetime, 
                        model, 
                        coords_top)
                    Fsums[j, θbin_index, ϕbin_index] += Fvalue
                end
            end
            return Fsums
        end
    end
    fetched_results = fetch.(tasks)
    # Initialize an array to hold the sum of `q` values in each bin
    Fsums_total = zeros(length(observation_energies), length(θbins)-1, length(ϕbins)-1)
    # Perform element-wise summation
    for result in fetched_results
        Fsums_total .+= result
    end
    for i in axes(Fsums_total, 2)
        Fsums_total[:, i, :] /= (cos(θbins[i]) - cos(θbins[i+1]))
    end
    D = configurations.max_radius
    Fsums_total /= 2π*D^2*period*(π/Nθ)*(2π/Nϕ)
    return midpoints(ϕbins), Fsums_total 
end

function observed_phase(pf::AbstractVector, period, ::SphericalTopology) 
    tf = pf[1]
    φf = pf[4]
    return mod(tf/period - φf/(2π), 1)    
end