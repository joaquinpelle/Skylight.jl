#Required
coordinates_topology(spacetime) = error("Coordinates topology not defined for this spacetime.")
metric!(g, position, spacetime, cache) = error("Metric not defined for this spacetime.")
metric!(g, position, spacetime::AbstractSpacetime, ::Nothing) = metric!(g, position, spacetime)

#Optional
allocate_cache(spacetime::AbstractSpacetime) = nothing
radius(position, spacetime) = error("Radius not defined for this spacetime.")
event_horizon_radius(spacetime) = error("Event horizon radius not defined for this spacetime.")
circular_geodesic_angular_speed(position, spacetime, rotation_sense) = error("Circular geodesic angular speed not defined for this spacetime.")
christoffel!(Γ, position, spacetime::AbstractSpacetime, ::Nothing) = christoffel!(Γ, position, spacetime)
allocate_christoffel_cache(spacetime::AbstractSpacetime) = AutoDiffChristoffelCache(spacetime)

"""
Computes the inverse of the given metric at the given position using a fast inversion
for 4x4 symmetric matrices.

Parameters:
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
axial_symmetry(spacetime::AbstractSpacetime) = isa(spherical_symmetry(spacetime), IsSphericallySymmetric) ? IsAxiallySymmetric() : IsNotAxiallySymmetric()

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
    cfg = ForwardDiff.JacobianConfig(spacetime_metric_field, Matrix{Float64}(undef, 4,4), Vector{Float64}(undef, 4), ForwardDiff.Chunk{4}())
    return AutoDiffChristoffelCache(spacetime_metric_field=spacetime_metric_field, spacetime_cache=spacetime_cache, cfg=cfg)
end

metric_field(spacetime::AbstractSpacetime, cache) = (g,position) -> metric!(g,position,spacetime,cache)

function metric(position, spacetime::AbstractSpacetime)
    g = zeros(4,4)
    cache = allocate_cache(spacetime)
    metric!(g, position, spacetime, cache)
    return g
end

function metric(position, spacetime::AbstractSpacetime, cache)
    g = zeros(4,4)
    metric!(g, position, spacetime, cache)
    return g
end

"""
Computes the volume element (square root of minus the determinant of the metric) at a given position 
using a fast determinant for 4x4 symmetric matrices.

Parameters:
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
    g = zeros(4,4)
    cache = allocate_cache(spacetime)
    metric!(g, position, spacetime, cache)
    return sqrt(-det4x4sym(g))
end

function metric_inverse(position, spacetime::AbstractSpacetime)
    ginv = zeros(4,4)
    g = zeros(4,4)
    cache = allocate_cache(spacetime)
    metric_inverse!(ginv, position, spacetime, g, cache)
    return ginv
end

function christoffel(position, spacetime::AbstractSpacetime)
    Γ = zeros(4,4,4)
    cache = allocate_christoffel_cache(spacetime)
    christoffel!(Γ, position, spacetime, cache)
    return Γ
end

function christoffel(position, spacetime::AbstractSpacetime, cache::AbstractChristoffelCache)
    Γ = zeros(4,4,4)
    christoffel!(Γ, position, spacetime, cache)
    return Γ
end

@inline sign(::ProgradeRotation) = 1.0
@inline sign(::RetrogradeRotation) = -1.0