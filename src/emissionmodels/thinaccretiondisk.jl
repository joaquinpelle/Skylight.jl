struct ThinAccretionDisk <: SurfaceEmissionModel end

function surface_function(position, model::ThinAccretionDisk, coord_system::CartesianKind)
    @views z = position[4]
    return z 
end

function surface_function(position, model::ThinAccretionDisk, coord_system::SphericalKind)
    @views θ = position[3]
    return θ-π/2 
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

function set_surface_differential!(covector, position, model::ThinAccretionDisk, coord_system::SphericalKind)

    @views begin
        t,r,θ,φ = position
    end

    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0

end    

