struct ThinAccretionDisk <: SurfaceEmissionModel end

function surface_function(position, model::ThinAccretionDisk, coord_system::CartesianKind)
    
    @views z = position[4]
    
    return z 

end

function set_surface_differential!(covector, position, model::ThinAccretionDisk, coord_system::CartesianKind)

    @views begin
        t,x,y,z = position
    end

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 0.0
    covector[4] = 1.0

end