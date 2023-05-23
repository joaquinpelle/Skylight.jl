abstract type AbstractRadiativeModel end
abstract type AbstractSurfaceEmissionModel <: AbstractRadiativeModel end

abstract type OpaqueInteriorSurfaceTrait end
struct IsOpaqueInteriorSurface <: OpaqueInteriorSurfaceTrait end
struct IsNotOpaqueInteriorSurface <: OpaqueInteriorSurfaceTrait end

#Set by default the radiative model to IsNotOpaqueInteriorSurface()
opaque_interior_surface_trait(::AbstractRadiativeModel) = IsNotOpaqueInteriorSurface()