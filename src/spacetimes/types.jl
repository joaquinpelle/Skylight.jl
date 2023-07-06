abstract type AbstractSpacetime end
abstract type AbstractBlackHoleSpacetime <: AbstractSpacetime end
abstract type AbstractRegularCompactObjectSpacetime <: AbstractSpacetime end

abstract type AbstractSpacetimeCache end
abstract type AbstractChristoffelCache end

abstract type AbstractCoordinatesTopology end
struct CartesianTopology <: AbstractCoordinatesTopology end
struct SphericalTopology <: AbstractCoordinatesTopology end

abstract type AbstractRotationSense end
struct ProgradeRotation <: AbstractRotationSense end
struct RetrogradeRotation <: AbstractRotationSense end

abstract type AbstractSpacetimeSymmetry end

abstract type Stationarity <: AbstractSpacetimeSymmetry end
struct IsStationary <: Stationarity end 
struct IsNotStationary <: Stationarity end 

abstract type SphericalSymmetry <: AbstractSpacetimeSymmetry end
struct IsSphericallySymmetric <: SphericalSymmetry end
struct IsNotSphericallySymmetric <: SphericalSymmetry end

abstract type AxialSymmetry <: AbstractSpacetimeSymmetry end
struct IsAxiallySymmetric <: AxialSymmetry end
struct IsNotAxiallySymmetric <: AxialSymmetry end

@with_kw mutable struct AutoDiffChristoffelCache{F,CA,CO} <: AbstractChristoffelCache
    g::Array{Float64,2} = zeros(4,4)
    ginv::Array{Float64,2} = zeros(4,4)
    ∂g::Array{Float64,3} = zeros(4,4,4)
    spacetime_metric_field::F
    spacetime_cache::CA
    cfg::CO
end