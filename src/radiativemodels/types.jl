abstract type AbstractRadiativeModel end
abstract type AbstractSurfaceEmissionModel <: AbstractRadiativeModel end

abstract type AbstractAccretionDisk <: AbstractSurfaceEmissionModel  end

abstract type OpaqueInteriorSurfaceTrait end
struct IsOpaqueInteriorSurface <: OpaqueInteriorSurfaceTrait end
struct IsNotOpaqueInteriorSurface <: OpaqueInteriorSurfaceTrait end

#Set by default the radiative model to IsNotOpaqueInteriorSurface()
opaque_interior_surface_trait(::AbstractRadiativeModel) = IsNotOpaqueInteriorSurface()

abstract type AbstractRadiativeProcess end
struct ThermalEmission <: AbstractRadiativeProcess end
struct Bremsstrahlung <: AbstractRadiativeProcess end

@with_kw struct Synchrotron{I} <: AbstractRadiativeProcess 
    T::Vector{Float64} = [5e8, 1e9, 2e9, 4e9, 8e9, 1e10, 3e10]
    α::I = LinearInterpolation([0.0431, 1.121, 1.180, 1.045, 0.9774, 0.9768, 0.9788], T)
    β::I = LinearInterpolation([10.44, -10.65, -4.008, -0.1897, 1.160, 1.095, 1.021], T)
    γ::I = LinearInterpolation([16.61, 9.169, 1.559, 0.0595, 0.2641, 0.8332, 1.031], T)
end