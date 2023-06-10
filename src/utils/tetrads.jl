function random_triad!(triad, time_vector, metric)
    @views begin
        triad_time_components = triad[1,:]
        triad_space_components = triad[2:4,:]  
    end
    fill!(triad_time_components,0.0)
    rand!(triad_space_components)
    orthonormalize!(triad, time_vector, metric)
    return nothing
end

function surface_adapted_triad!(triad, time_vector, metric, metric_inverse, position, model::AbstractRadiativeModel, coords_top)
    @views begin
        triad_time_components = triad[1,:]
        normal = triad[:,1]
        dyad = triad[:,2:3]
        dyad_space_components = dyad[2:4,:]
    end
    
    fill!(triad_time_components,0.0)
    rand!(dyad_space_components)
    unit_surface_normal!(normal, position, metric, metric_inverse, model, coords_top) 
    orthonormalize!(triad, time_vector, metric)
    return nothing
end

function spherical_like_triad!(triad, position, time_vector, metric, ::CartesianTopology)
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    
    fill!(triad, 0.0)    
    @views begin
        e1 = triad[2:4,1]
        e2 = triad[2:4,2]
        e3 = triad[2:4,3]
    end

    e1[1] = -x/r
    e1[2] = -y/r
    e1[3] = -z/r

    if !(x==y==0.0)
        e2[1] = -y/r
        e2[2] =  x/r
    else
        e2[1] = 0.0
        e2[2] = 1.0
    end

    e3 .= -cross(e1,e2)
    orthonormalize!(triad, time_vector, metric)
end

function spherical_like_triad!(triad, position, time_vector, metric, ::SphericalTopology)
    fill!(triad, 0.0)    
    @views begin
        e1 = triad[2:4,1]
        e2 = triad[2:4,2]
        e3 = triad[2:4,3]
    end
    e1[1] = -1.0
    e2[3] =  1.0
    e3[2] = -1.0
    orthonormalize!(triad, time_vector, metric)
end

"""
    orthonormalize!(dyad::Matrix, time_vector::Vector, space_vector::Vector, metric)

    Orthonormalize the dyad `dyad` with respect to the metric `metric` and the vectors `time_vector` and `space_vector`.
    Both `time_vector` and `space_vector` are assumed to be normalized and orthogonal to each other.
"""
function orthonormalize!(dyad, time_vector, space_vector, metric)
    @views begin
        e1 = dyad[:,1]
        e2 = dyad[:,2]
    end

    e1 .+= time_vector*scalar_product(e1,time_vector,metric) - space_vector*scalar_product(e1,space_vector,metric)
    normalize_spacelike!(e1,metric)

    e2 .+= time_vector*scalar_product(e2,time_vector,metric) - space_vector*scalar_product(e2,space_vector,metric) - e1*scalar_product(e2,e1,metric)
    normalize_spacelike!(e2,metric)
    return nothing
end

"""
    orthonormalize!(triad::Matrix, time_vector::Vector, metric)

    Orthonormalize the triad `triad` with respect to the metric `metric` and the vector `time_vector`.
    `time_vector` is assumed to be normalized.
"""
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
    return nothing
end

function unit_time_components!(kμ)
    kμ[1,:] .= 1.0
    return nothing
end

function negative_unit_time_components!(kμ)
    kμ[1,:] .= -1.0
    return nothing
end

function coordinate_components_from_tetrad_components!(vectors, tetrad)
    for i in axes(vectors,2)
        vectors[:,i] = tetrad*vectors[:,i]  
    end
    return nothing
end