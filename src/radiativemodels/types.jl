abstract type AbstractRadiativeModel end

abstract type SurfaceEmissionModel <: AbstractRadiativeModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end
