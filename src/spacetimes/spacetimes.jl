#Required
function coordinates_topology(spacetime)
    error("Coordinates topology not defined for this spacetime.")
end
metric!(g, position, spacetime, cache) = error("Metric not defined for this spacetime.")
function metric!(g, position, spacetime::AbstractSpacetime, ::Nothing)
    metric!(g, position, spacetime)
end

#Optional
allocate_cache(::AbstractSpacetime) = nothing
radius(position, spacetime) = error("Radius not defined for this spacetime.")
function event_horizon_radius(spacetime)
    error("Event horizon radius not defined for this spacetime.")
end
function circular_geodesic_angular_speed(position, spacetime, rotation_sense)
    error("Circular geodesic angular speed not defined for this spacetime.")
end
function christoffel!(Γ, position, spacetime::AbstractSpacetime, ::Nothing)
    christoffel!(Γ, position, spacetime)
end
function allocate_christoffel_cache(spacetime::AbstractSpacetime)
    AutoDiffChristoffelCache(spacetime)
end

"""
Computes the inverse of the given metric at the given position using a fast inversion
for 4x4 symmetric matrices.

Arguments:
- ginv: mutable array of size (4,4) to store the resulting inverse metric.
- position: tuple of four numbers representing a position in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) to store the metric evaluated at the given position.

Returns: nothing.
"""
function metric_inverse!(ginv, position, spacetime::AbstractSpacetime, g, cache)
    metric!(g, position, spacetime, cache)
    inv4x4sym!(ginv, g)
    return nothing
end

#By default we set non-stationarity and non-spherical symmetry  
stationarity(::AbstractSpacetime) = IsNotStationary()
spherical_symmetry(::AbstractSpacetime) = IsNotSphericallySymmetric()
#For z symmetry we check spherical symmetry by default first. Thus if a spacetime is declared spherically symmetric it's automatically axially symmetric.
function axial_symmetry(spacetime::AbstractSpacetime)
    isa(spherical_symmetry(spacetime), IsSphericallySymmetric) ? IsAxiallySymmetric() :
    IsNotAxiallySymmetric()
end

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

function metric_field(spacetime::AbstractSpacetime, cache)
    (g, position) -> metric!(g, position, spacetime, cache)
end

function metric(position, spacetime::AbstractSpacetime)
    g = zeros(4, 4)
    cache = allocate_cache(spacetime)
    metric!(g, position, spacetime, cache)
    return g
end

function metric(position, spacetime::AbstractSpacetime, cache)
    g = zeros(4, 4)
    metric!(g, position, spacetime, cache)
    return g
end

"""
Computes the volume element (square root of minus the determinant of the metric) at a given position 
using a fast determinant for 4x4 symmetric matrices.

Arguments:
- position: tuple of four numbers representing a position in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) to store the metric evaluated at the given position.

Returns: the volume element.
"""
function volume_element(position, spacetime::AbstractSpacetime, g, cache)
    metric!(g, position, spacetime, cache)
    return sqrt(-det4x4sym(g))
end

function volume_element(position, spacetime::AbstractSpacetime)
    g = zeros(4, 4)
    cache = allocate_cache(spacetime)
    metric!(g, position, spacetime, cache)
    return sqrt(-det4x4sym(g))
end

function metric_inverse(position, spacetime::AbstractSpacetime)
    ginv = zeros(4, 4)
    g = zeros(4, 4)
    cache = allocate_cache(spacetime)
    metric_inverse!(ginv, position, spacetime, g, cache)
    return ginv
end

function christoffel(position, spacetime::AbstractSpacetime)
    Γ = zeros(4, 4, 4)
    cache = allocate_christoffel_cache(spacetime)
    christoffel!(Γ, position, spacetime, cache)
    return Γ
end

function christoffel(position,
    spacetime::AbstractSpacetime,
    cache::AbstractChristoffelCache)
    Γ = zeros(4, 4, 4)
    christoffel!(Γ, position, spacetime, cache)
    return Γ
end

@inline sign(::ProgradeRotation) = 1.0
@inline sign(::RetrogradeRotation) = -1.0

""" 
    equatorial_ring_areas(edges, spacetime)

Approximate areas of the equatorial rings delimited by `edges` in `spacetime`. The spacetime
must be axisymmetric. The areas are computed using the formula
`2π*sqrt(g[2,2]*g[4,4])*Δr`, where `g` is the metric evaluated at the center of
the ring.
"""
function equatorial_ring_areas(edges::AbstractVector, spacetime::AbstractSpacetime)
    is_axisymmetric(spacetime) || throw(ArgumentError("Spacetime must be axisymmetric"))
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