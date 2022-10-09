function set_local_triad!(triad, position, time_vector, metric, metric_inverse, model, coord_system)

    set_random_triad!(triad) 
    orthonormalize!(triad, time_vector, metric)

end

function set_local_triad!(triad, position, time_vector, metric, metric_inverse, model::SurfaceEmissionModel, coord_system)

    @views begin
        e1 = triad[:,1]
        dyad = triad[:,2:3]
    end
    
    set_unit_surface_normal!(e1, position, metric, metric_inverse, model, coord_system) 
    set_random_tangent_dyad!(dyad)
    orthonormalize!(triad, time_vector, metric)

end

function set_random_tangent_dyad!(dyad)

    fill!(dyad,0.0)
    dyad[2:4,:] = rand(3,2)

end

function set_random_triad!(triad)

    fill!(triad,0.0)
    triad[2:4,:] = rand(3,3)

end

function orthonormalize!(triad, time_vector, metric)

    @views begin
        e1 = triad[:,1]
        e2 = triad[:,2]
        e3 = triad[:,3]
    end

    e1 .+= time_vector*scalar_product(e1,time_vector,metric)
    normalize_spacelike!(e1,metric)

    e2 .+= time_vector*scalar_product(e2,time_vector,metric) - e1*scalar_product(e2,e1,metric)
    normalize_spacelike!(e2,metric)

    e3 .+= time_vector*scalar_product(e3,time_vector,metric) - e1*scalar_product(e3,e1,metric) - e2*scalar_product(e3,e2,metric)
    normalize_spacelike!(e3,metric)

end

function set_coordinate_components_from_tetrad_components!(kμ, tetrad)

    nvectors = size(kμ,2)

    for i in 1:nvectors
        
        #kμ[:,i] = tetrad[:,1]*kμ[1,i] + tetrad[:,2]*kμ[2,i] + tetrad[:,3]*kμ[3,i] + tetrad[:,4]*kμ[4,i]  
        kμ[:,i] = tetrad*kμ[:,i]  
    
    end

end