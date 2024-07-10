"""
    AbstractSpacetime
"""
abstract type AbstractSpacetime end

"""
    AbstractBlackHoleSpacetime <: AbstractSpacetime
"""
abstract type AbstractBlackHoleSpacetime <: AbstractSpacetime end

"""
    AbstractRegularCompactObjectSpacetime <: AbstractSpacetime
"""
abstract type AbstractRegularCompactObjectSpacetime <: AbstractSpacetime end

"""
AbstractSpacetimeCache
"""
abstract type AbstractSpacetimeCache end


"""
AbstractSpacetimeCache
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
    spacetime_metric_closure::F
    spacetime_cache::CA
    cfg::CO
end

"""
    AbstractCoordinatesTopology
"""
abstract type AbstractCoordinatesTopology end

"""
    CartesianTopology <: AbstractCoordinatesTopology
"""
struct CartesianTopology <: AbstractCoordinatesTopology end

"""
    SphericalTopology <: AbstractCoordinatesTopology
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

abstract type AbstractSymmetry end

abstract type Stationarity <: AbstractSymmetry end
struct IsStationary <: Stationarity end
struct IsNotStationary <: Stationarity end

abstract type SphericalSymmetry <: AbstractSymmetry end
struct IsSphericallySymmetric <: SphericalSymmetry end
struct IsNotSphericallySymmetric <: SphericalSymmetry end

abstract type AxialSymmetry <: AbstractSymmetry end
struct IsAxisymmetric <: AxialSymmetry end
struct IsNotAxisymmetric <: AxialSymmetry end

abstract type HelicalSymmetry <: AbstractSymmetry end
struct IsHelicallySymmetric <: HelicalSymmetry end
struct IsNotHelicallySymmetric <: HelicalSymmetry end