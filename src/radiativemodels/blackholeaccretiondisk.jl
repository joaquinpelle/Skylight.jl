@with_kw struct BlackHoleAccretionDisk <: AbstractSurfaceEmissionModel 

    inner_radius::Float64
    outer_radius::Float64 

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

function is_final_position_at_source(position, spacetime, model::BlackHoleAccretionDisk)

    r = position[2]

    return (r >= model.inner_radius) && (r <= model.outer_radius)

end