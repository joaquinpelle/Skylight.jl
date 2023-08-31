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

Abstract type for cache objects to be used as temporary storage in Christoffel symbol calculations.
"""
abstract type AbstractChristoffelCache end

"""
    AutoDiffChristoffelCache{F, CA, CO} <: AbstractChristoffelCache

Cache object for temporary storage in Christoffel symbol calculation via automatic differentiation. 

# Constructor

```
AutoDiffChristoffelCache(spacetime::AbstractSpacetime)
```
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

"""
    AbstractRotationSense

Abstract type for representing the rotation sense of a circular orbit in a spacetime.
"""
abstract type AbstractRotationSense end

"""
    ProgradeRotation <: AbstractRotationSense

Prograde rotation sense for a circular orbit in a spacetime.

# Constructor

```
ProgradeRotation()
```
"""
struct ProgradeRotation <: AbstractRotationSense end

"""
    RetrogradeRotation <: AbstractRotationSense

Retrograde rotation sense for a circular orbit in a spacetime.

# Constructor

```
RetrogradeRotation()
```
"""
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