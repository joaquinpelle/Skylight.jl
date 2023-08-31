"""
    AbstractSpacetime

Supertype for representing the geometrical structure of spacetime. Specific implementations of spacetime, such as black holes and regular compact objects, are subtypes of this abstract type.
"""
abstract type AbstractSpacetime end

"""
    AbstractBlackHoleSpacetime <: AbstractSpacetime

Supertype for representing spacetimes containing black holes. 
"""
abstract type AbstractBlackHoleSpacetime <: AbstractSpacetime end

"""
    AbstractRegularCompactObjectSpacetime <: AbstractSpacetime

Supertype for representing spacetimes containing compact objects without horizons like boson stars. 
"""
abstract type AbstractRegularCompactObjectSpacetime <: AbstractSpacetime end

"""
AbstractSpacetimeCache

Abstract type for caching spacetime-related computations. This can be used as scratch memory in calculations involving the spacetime.
"""
abstract type AbstractSpacetimeCache end


"""
AbstractSpacetimeCache

Abstract type for caching Christoffel symbol-related calculations. This can be used as scratch memory in calculations in the `christoffel!` function.
"""
abstract type AbstractChristoffelCache end

"""
AutoDiffChristoffelCache{F, CA, CO} <: AbstractChristoffelCache

Mutable structure for caching values needed when computing Christoffel symbols using automatic differentiation. It includes the metric tensor `g`, its inverse `ginv`, derivatives of the metric `∂g`, a spacetime metric field `spacetime_metric_field`, a cache for spacetime calculations `spacetime_cache`, and a Jacobian configuration `cfg`.

# Fields
- `g::Array{Float64, 2}`: The metric tensor (default: zeros(4, 4)).
- `ginv::Array{Float64, 2}`: The inverse of the metric tensor (default: zeros(4, 4)).
- `∂g::Array{Float64, 3}`: Derivatives of the metric tensor (default: zeros(4, 4, 4)).
- `spacetime_metric_field::F`: A field representing the spacetime metric.
- `spacetime_cache::CA`: A cache for spacetime-related computations.
- `cfg::CO`: Configuration object for additional settings.
"""
@with_kw mutable struct AutoDiffChristoffelCache{F, CA, CO} <: AbstractChristoffelCache
    g::Array{Float64, 2} = zeros(4, 4)
    ginv::Array{Float64, 2} = zeros(4, 4)
    ∂g::Array{Float64, 3} = zeros(4, 4, 4)
    spacetime_metric_field::F
    spacetime_cache::CA
    cfg::CO
end

"""
    AbstractCoordinatesTopology

Abstract type for representing the topology of the coordinates of a spacetime.
"""
abstract type AbstractCoordinatesTopology end

"""
    CartesianTopology <: AbstractCoordinatesTopology

Cartesian coordinates topology.
"""
struct CartesianTopology <: AbstractCoordinatesTopology end

"""
    SphericalTopology <: AbstractCoordinatesTopology

Spherical coordinates topology.
"""
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