@with_kw struct BlackHoleAccretionDisk{T} <: AbstractSurfaceEmissionModel 

    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation() 

end

function set_surface_differential!(covector, position, ::BlackHoleAccretionDisk, ::CartesianTopology)

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 0.0
    covector[4] = 1.0

end

function set_surface_differential!(covector, position, ::BlackHoleAccretionDisk, ::SphericalTopology)

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0

end    

function set_emitter_four_velocity!(vector, position, metric, spacetime::AbstractKerrSpacetime, model::BlackHoleAccretionDisk, coords_top)
    angular_speed = circular_geodesic_angular_speed(position, spacetime, model.rotation_sense)
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coords_top)
end

function is_final_position_at_source(position, spacetime::AbstractKerrSpacetime, model::BlackHoleAccretionDisk)
    r = kerr_radius(position, spacetime)
    return (r >= model.inner_radius) && (r <= model.outer_radius)
end