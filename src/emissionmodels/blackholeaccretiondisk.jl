export NovikovThorneDisk

@with_kw struct NovikovThorneDisk <: BlackHoleAccretionDisk

    inner_radius::Float64
    outer_radius::Float64 
    rbound::Float64 = 1e-6

end

function set_surface_differential!(covector, position, model::BlackHoleAccretionDisk, coord_system::CartesianClass)

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 0.0
    covector[4] = 1.0

end

function set_surface_differential!(covector, position, model::BlackHoleAccretionDisk, coord_system::SphericalClass)

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0

end    

