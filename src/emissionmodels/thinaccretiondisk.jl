struct ThinAccretionDisk <: SurfaceEmissionModel end

function set_surface_differential!(covector, position, model::ThinAccretionDisk, coord_system::CartesianKind)

    @views begin
        t,x,y,z = position
    end

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 0.0
    covector[4] = 1.0

end

function set_surface_differential!(covector, position, model::ThinAccretionDisk, coord_system::SphericalKind)

    @views begin
        t,r,θ,φ = position
    end

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0

end    

