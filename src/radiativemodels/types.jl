abstract type RadiativeModel end

abstract type SurfaceEmissionModel <: RadiativeModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end
