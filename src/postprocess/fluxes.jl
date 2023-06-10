function fluxes(Iobs, configurations, ::ImagePlane)
    Fobs = copy(Iobs)
    normalize_by_pixel_area!(Fobs, configurations)
    normalize_by_image_plane_distance!(Fobs, configurations)
end

"""For monochromatic / bolometric fluxes"""
function fluxes(I::AbstractVector, configurations, camera::PinholeCamera, initial_data::AbstractMatrix, output_data::AbstractMatrix; observer_four_velocity=nothing, flux_direction=nothing)

    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    size(initial_data, 2) == length(I) || throw(DimensionMismatch("The number of rays in the initial and output data must be the same as the length of I."))
    number_of_pixels(camera) == length(I) || throw(DimensionMismatch("The number of pixels in the camera must be the same as the length of I."))

    spacetime = configurations.spacetime
    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    flux_direction!(cache, flux_direction, camera, spacetime) 

    dΩ = pixel_solid_angles(camera)
    F = zero(I)
    @inbounds begin
        for i in eachindex(I)

            @views begin 
                ki = initial_data[5:8,i]
                pf = output_data[1:4,i]
            end

            if !is_final_position_at_source(pf, spacetime, model)
                continue
            end
            nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
            nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
            F[i] = nu*nn*I[i]*dΩ[i]        
        end
    end
    return F
end

"""For multi-energies fluxes"""
function fluxes(I::AbstractMatrix, camera::PinholeCamera, initial_data::AbstractMatrix, output_data::AbstractMatrix, spacetime::AbstractSpacetime; observer_four_velocity=nothing, flux_direction=nothing)

    same_size(initial_data, output_data) || throw(DimensionMismatch("initial_data and output_data must have the same size."))
    size(initial_data, 2) == length(I) || throw(DimensionMismatch("The number of rays in the initial and output data must be the same as the length of I."))
    number_of_pixels(camera) == length(I) || throw(DimensionMismatch("The number of pixels in the camera must be the same as the length of I."))

    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    flux_direction!(cache, flux_direction, camera, spacetime) 

    dΩ = pixel_solid_angles(camera)
    Fobs = zeros(Iobs)
    @inbounds begin
        for i in eachindex(I)

            @views begin 
                ki = initial_data[5:8,i]
                pf = output_data[1:4,i]
            end

            if !is_final_position_at_source(pf, spacetime, model)
                continue
            end

            nu = scalar_product(ki, cache.observer_four_velocity, cache.observer_metric)
            nn = scalar_product(ki, cache.flux_direction, cache.observer_metric)
            
            Fobs[:,i] = nu*nn*I[:,i]*dΩ[i]        
        end
    end
    return Fobs
end