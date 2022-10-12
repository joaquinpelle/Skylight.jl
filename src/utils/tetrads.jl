function set_triad!(triad, time_vector, metric)

    @views begin
        triad_time_components = triad[1,:]
        triad_space_components = triad[2:4,:]  
    end

    fill!(triad_time_components,0.0)
    rand!(triad_space_components)
    orthonormalize!(triad, time_vector, metric)

end

function set_surface_adapted_triad!(triad, time_vector, position, metric, metric_inverse, model, coord_system)

    @views begin
        triad_time_components = triad[1,:]
        normal = triad[:,1]
        dyad = triad[:,2:3]
        dyad_space_components = dyad[2:4,:]
    end
    
    fill!(triad_time_components,0.0)
    rand!(dyad_space_components)

    set_unit_surface_normal!(normal, position, metric, metric_inverse, model, coord_system) 
    orthonormalize!(triad, time_vector, metric)

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

function set_coordinate_components_from_tetrad_components!(vectors, tetrad)

    nvectors = size(vectors,2)

    for i in 1:nvectors
        
        vectors[:,i] = tetrad*vectors[:,i]  
    
    end

end