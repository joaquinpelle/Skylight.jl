#Required interface
function coordinates_topology(spacetime)
    error("Coordinates topology not defined for this spacetime.")
end

"""
    metric!(g::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractSpacetimeCache})

Evaluate the spacetime metric at the given position and store the result in `g` using `cache` for temporary storage.  
"""
metric!(g, position, spacetime, cache) = error("Metric not defined for this spacetime.")

"""
    allocate_cache(spacetime::AbstractSpacetime)

Allocate a cache object for the given spacetime. The cache object is used to store temporary data in spacetime-related calculations.
"""
allocate_cache(::AbstractSpacetime) = nothing

include("coordinatealias.jl")
include("autodiff.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("frkerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("rar.jl")
include("numerical.jl")
include("constantsmotion.jl")

"""
    metric(position::AbstractVector, spacetime::AbstractSpacetime)

Evaluate the spacetime metric at the given position and return the result.
"""
function metric(position::AbstractVector, spacetime::AbstractSpacetime)
    g = zeros(4, 4)
    cache = allocate_cache(spacetime)
    metric!(g, position, spacetime, cache)
    return g
end

"""
    metric(position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractSpacetimeCache})

Evaluate the spacetime metric at the given position and return the result using a cache object as temporary storage.
The cache must be pre-allocated as `cache = allocate_cache(spacetime)`. 

See also [`allocate_cache`](@ref). 
"""
function metric(position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractSpacetimeCache})
    g = zeros(4, 4)
    metric!(g, position, spacetime, cache)
    return g
end

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, ::Nothing)
    metric!(g, position, spacetime)
end

function metric_field(spacetime::AbstractSpacetime)
    (g, position) -> metric!(g, position, spacetime)
end

function metric_field(spacetime::AbstractSpacetime, cache::AbstractSpacetimeCache)
    (g, position) -> metric!(g, position, spacetime, cache)
end

function metric_field(spacetime::AbstractSpacetime, ::Nothing)
    metric_field(spacetime)
end

"""
    metric_inverse(position::AbstractVector, spacetime::AbstractSpacetime)

Evaluate the inverse of the metric at the given position and return the result.
"""
function metric_inverse(position::AbstractVector, spacetime::AbstractSpacetime)
    ginv = zeros(4, 4)
    g = zeros(4, 4)
    cache = allocate_cache(spacetime)
    metric_inverse!(ginv, position, spacetime, g, cache)
    return ginv
end

"""
    metric_inverse!(ginv::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})

Evaluate the inverse of the metric at the given position  and store the result in the given array, using `g` and `cache` for temporary storage, which
must be preallocated as `g = zeros(4,4)` and `cache = allocate_cache(spacetime)`.

See also [`allocate_cache`](@ref). 
"""
function metric_inverse!(ginv::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
    metric!(g, position, spacetime, cache)
    inv4x4sym!(ginv, g)
    return nothing
end

function metric_inverse!(ginv::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, ::Nothing)
    metric_inverse!(ginv, position, spacetime, g)
end

function metric_inverse!(ginv::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix)
    metric!(g, position, spacetime)
    inv4x4sym!(ginv, g)
    return nothing
end

@doc raw"""
    volume_element(position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})

Compute $\sqrt{-\text{det}(g)}$ at the given position, using `g` and `cache` for temporary storage, which must be preallocated as `g = zeros(4,4)` and `cache = allocate_cache(spacetime)`. 

See also [`allocate_cache`](@ref). 
"""
function volume_element(position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
    metric!(g, position, spacetime, cache)
    return sqrt(-det4x4sym(g))
end

@doc raw"""
    volume_element(position::AbstractVector, spacetime::AbstractSpacetime)

Compute $\sqrt{-\text{det}(g)}$ at the given position 
"""
function volume_element(position::AbstractVector, spacetime::AbstractSpacetime)
    g = zeros(4, 4)
    cache = allocate_cache(spacetime)
    metric!(g, position, spacetime, cache)
    return sqrt(-det4x4sym(g))
end

@doc raw"""
    christoffel(position::AbstractVector, spacetime::AbstractSpacetime)

Evaluate the Christoffel symbols of the second kind at the given position and return the result as a `4x4x4` array 
`Γ`, where `Γ[α,μ,ν]` corresponds to 

``\Gamma^\alpha_{\mu \nu} = \frac{1}{2} g^{\alpha \rho}(\partial_\mu g_{\rho \nu} + \partial_\nu g_{\rho \mu} - \partial_\rho g_{\mu \nu}).``
"""
function christoffel(position::AbstractVector, spacetime::AbstractSpacetime)
    Γ = zeros(4, 4, 4)
    cache = allocate_christoffel_cache(spacetime)
    christoffel!(Γ, position, spacetime, cache)
    return Γ
end

function christoffel(position::AbstractVector,
    spacetime::AbstractSpacetime,
    cache::Union{Nothing,AbstractChristoffelCache})
    Γ = zeros(4, 4, 4)
    christoffel!(Γ, position, spacetime, cache)
    return Γ
end

@doc raw"""
    christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing,AbstractChristoffelCache})

Evaluate the Christoffel symbols of the second kind at the given position and store the result in
`Γ`, where `Γ[α,μ,ν]` corresponds to 

``\Gamma^\alpha_{\mu \nu} = \frac{1}{2} g^{\alpha \rho}(\partial_\mu g_{\rho \nu} + \partial_\nu g_{\rho \mu} - \partial_\rho g_{\mu \nu}).``

The cache must be preallocated as `cache = allocate_christoffel_cache(spacetime)`.

See also [`allocate_christoffel_cache`](@ref).
"""
function christoffel!(::AbstractArray, ::AbstractVector, ::AbstractSpacetime, ::Union{Nothing,AbstractChristoffelCache})
    error("Christoffel symbols not defined for this spacetime. This is weird because the behaviour should default to automatic differentiation.
    Contact the author for help.")
end

function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::AbstractSpacetime, ::Nothing)
    christoffel!(Γ, position, spacetime)
end

"""
    allocate_christoffel_cache(spacetime::AbstractSpacetime)

Allocate a cache object for Christoffel symbols calculations.

See also [`christoffel!`](@ref). 
""" 
function allocate_christoffel_cache(spacetime::AbstractSpacetime)
    AutoDiffChristoffelCache(spacetime)
end

function AutoDiffChristoffelCache(spacetime::AbstractSpacetime)
    spacetime_cache = allocate_cache(spacetime)
    spacetime_metric_field = metric_field(spacetime, spacetime_cache)
    cfg = ForwardDiff.JacobianConfig(spacetime_metric_field,
        Matrix{Float64}(undef, 4, 4),
        Vector{Float64}(undef, 4),
        ForwardDiff.Chunk{4}())
    return AutoDiffChristoffelCache(spacetime_metric_field = spacetime_metric_field,
        spacetime_cache = spacetime_cache,
        cfg = cfg)
end

"""
    radius(position, spacetime)

Radius of the given position in the given spacetime. 
"""
radius(position, spacetime) = error("Radius not defined for this spacetime.")

"""
    event_horizon_radius(spacetime)

Radius of the event horizon of the given spacetime. Defined only for black hole spacetimes.
"""
function event_horizon_radius(spacetime)
    error("Event horizon radius not defined for this spacetime.")
end

"""
    isco_radius(spacetime, rotation_sense::AbstractRotationSense)

Radius of the innermost stable circular orbit in the given spacetime and rotation sense (prograde or 
retrograde), if defined.
"""
function isco_radius(spacetime, rotation_sense)
    error("Isco radius radius not defined for this spacetime.")
end

"""
    mbco_radius(spacetime, rotation_sense::AbstractRotationSense)

Radius of the marginally bound circular orbit in the given spacetime and rotation sense (prograde or
retrograde), if defined.
"""
function mbco_radius(spacetime, rotation_sense)
    error("Marginally bound circular orbit not defined for this spacetime.")
end

"""
    circular_geodesic_angular_speed(position, spacetime, rotation_sense::AbstractRotationSense)

Angular speed of a circular geodesic at the given position in the given spacetime for a particle rotating
in the given sense.
""" 
function circular_geodesic_angular_speed(position, spacetime, rotation_sense)
    error("Circular geodesic angular speed not defined for this spacetime.")
end

#Other functions
#By default we set non-stationarity and non-spherical symmetry  
"""
    stationarity(spacetime)

Return `IsStationary()` if the spacetime is stationary, `IsNotStationary()` otherwise.
"""
stationarity(::AbstractSpacetime) = IsNotStationary()

"""
    spherical_symmetry(spacetime)

Return `IsSphericallySymmetric()` if the spacetime is spherically symmetric, `IsNotSphericallySymmetric()` otherwise.
"""
spherical_symmetry(::AbstractSpacetime) = IsNotSphericallySymmetric()
#For z symmetry we check spherical symmetry by default first. Thus if a spacetime is declared spherically symmetric it's automatically axially symmetric.
"""
    axial_symmetry(spacetime)

Return `IsAxiallySymmetric()` if the spacetime is axially symmetric, `IsNotAxiallySymmetric()` otherwise.
"""
function axial_symmetry(spacetime::AbstractSpacetime)
    isa(spherical_symmetry(spacetime), IsSphericallySymmetric) ? IsAxiallySymmetric() :
    IsNotAxiallySymmetric()
end

@inline sign(::ProgradeRotation) = 1.0
@inline sign(::RetrogradeRotation) = -1.0

"""
    equatorial_ring_areas(edges::AbstractVector, spacetime::AbstractSpacetime)

 
    equatorial_ring_areas(edges, spacetime)

Approximate areas of the equatorial rings delimited by `edges` in `spacetime`. The spacetime
must be stationary and axisymmetric. The areas are computed using the formula
`2π*sqrt(g[2,2]*g[4,4])*Δr`, where `g` is the metric evaluated at the center of
the ring.
"""
function equatorial_ring_areas(edges::AbstractVector, spacetime::AbstractSpacetime)
    is_stationary(spacetime) || throw(ArgumentError("Spacetime must be axisymmetric"))
    is_axisymmetric(spacetime) || throw(ArgumentError("Spacetime must be stationary"))
    coords_top = coordinates_topology(spacetime)
    centers = midpoints(edges)
    Δr = widths(edges)
    areas = zeros(length(centers))
    g = zeros(4,4)
    position = zeros(4)
    for (i,r) in enumerate(radii)
        position = equatorial_position!(position, r, coords_top)
        metric!(g, position, spacetime)
        areas[i] = 2π*sqrt(g[2,2]*g[4,4])*Δr[i]
    end
    return areas
end