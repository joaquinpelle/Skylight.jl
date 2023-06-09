function fluxes(Iobs, camera::ImagePlane)
    Fobs = copy(Iobs)
    normalize_by_pixel_area!(Fobs, camera)
    normalize_by_image_plane_distance!(Fobs, camera)
end

"""For monochromatic / bolometric fluxes"""
function fluxes(I::AbstractVector, camera::PinholeCamera, initial_data::AbstractMatrix, output_data::AbstractMatrix, spacetime::AbstractSpacetime; observer_four_velocity=nothing, flux_direction=nothing)

    Nrays = number_of_initial_conditions(camera)
    @assert size(I,1)==Nrays "The number of rays in the initial data and the number of rays in the image plane must be the same."

    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    flux_direction!(cache, flux_direction, camera, spacetime) 

    dΩ = pixel_solid_angles(camera)
    F = zeros(I)
    @inbounds begin
        for i in 1:Nrays

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

    Nrays = number_of_initial_conditions(camera)
    @assert size(I,2)==Nrays "The number of rays in the initial data and the number of rays in the image plane must be the same."

    cache = postprocess_cache(camera)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    flux_direction!(cache, flux_direction, camera, spacetime) 

    dΩ = pixel_solid_angles(camera)
    Fobs = zeros(Iobs)
    @inbounds begin
        for i in 1:Nrays

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