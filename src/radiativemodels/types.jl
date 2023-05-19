abstract type AbstractRadiativeModel end

abstract type AbstractSurfaceEmissionModel <: AbstractRadiativeModel end

abstract type AbstractOpaqueInteriorSurfaceEmissionModel <: AbstractSurfaceEmissionModel end
