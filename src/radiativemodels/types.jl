abstract type AbstractRadiativeProcess end
struct Bremsstrahlung <: AbstractRadiativeProcess end
struct Synchrotron <: AbstractRadiativeProcess end
struct InverseCompton <: AbstractRadiativeProcess end
struct ThermalEmission <: AbstractRadiativeProcess end

abstract type AbstractRadiativeModel end
abstract type AbstractSurfaceEmissionModel <: AbstractRadiativeModel end

abstract type AbstractAccretionDisk <: AbstractSurfaceEmissionModel  end

abstract type OpaqueInteriorSurfaceTrait end
struct IsOpaqueInteriorSurface <: OpaqueInteriorSurfaceTrait end
struct IsNotOpaqueInteriorSurface <: OpaqueInteriorSurfaceTrait end

#Set by default the radiative model to IsNotOpaqueInteriorSurface()
opaque_interior_surface_trait(::AbstractRadiativeModel) = IsNotOpaqueInteriorSurface()