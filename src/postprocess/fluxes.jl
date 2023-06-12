function fluxes!(F, configurations, ::ImagePlane)
    normalize_by_pixel_area!(F, configurations)
    normalize_by_image_plane_distance!(F, configurations)
    return F
end

"""For monochromatic / bolometric fluxes"""
function fluxes!(F::AbstractVector, configurations, camera::PinholeCamera, initial_data::AbstractMatrix, output_data::AbstractMatrix; observer_four_velocity=nothing, flux_direction=nothing)

    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    size(initial_data, 2) == size(I,1) || throw(DimensionMismatch("The number of rays in the initial and output data must be the same as the length of I."))
    number_of_pixels(camera) == size(I,1) || throw(DimensionMismatch("The number of pixels in the camera must be the same as the length of I."))

    spacetime = configurations.spacetime
    threads_cache = postprocess_cache(camera)

    for cache in threads_cache
        observer_metric!(cache, camera.position, spacetime)
        observer_four_velocity!(cache, observer_four_velocity) 
        flux_direction!(cache, flux_direction, camera, spacetime) 
    end

    dΩ = pixel_solid_angles(camera)
    @inbounds begin
        @threads for i in axes(initial_data,2)
            cache = threads_cache[threadid()]

            @views begin 
                ki = initial_data[5:8,i]
                pf = output_data[1:4,i]
            end
            if !is_final_position_at_source(pf, spacetime, model)
                continue
            end
            nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
            nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
            F[i] *= nu*nn*dΩ[i]        
        end
    end
    return F
end

"""For multi-energies fluxes"""
function fluxes!(F::AbstractMatrix, configurations, camera::PinholeCamera, initial_data::AbstractMatrix, output_data::AbstractMatrix; observer_four_velocity=nothing, flux_direction=nothing)

    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    sane_size(initial_data, I, 2) || throw(DimensionMismatch("The number of rays in the initial and output data must be the same as the length of I."))
    number_of_pixels(camera) == size(I,2) || throw(DimensionMismatch("The number of pixels in the camera must be the same as the length of I."))

    spacetime = configurations.spacetime
    threads_cache = postprocess_cache(camera)

    for cache in threads_cache
        observer_metric!(cache, camera.position, spacetime)
        observer_four_velocity!(cache, observer_four_velocity) 
        flux_direction!(cache, flux_direction, camera, spacetime) 
    end

    dΩ = pixel_solid_angles(camera)
    @inbounds begin
        @threads for i in axes(initial_data,2)
            cache = threads_cache[threadid()]
            @views begin 
                ki = initial_data[5:8,i]
                pf = output_data[1:4,i]
            end

            if !is_final_position_at_source(pf, spacetime, model)
                continue
            end
            nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
            nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
            F[:,i] *= nu*nn*dΩ[i]        
        end
    end
    return Fobs
end