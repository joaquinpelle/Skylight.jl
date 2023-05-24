#Required
coordinates_topology(::AbstractSpacetime) = error("Coordinates topology not defined for this spacetime.")
set_metric!(::AbstractSpacetime) = error("Metric not defined for this spacetime.")
set_christoffel!(Γ, position, ::AbstractSpacetime, cache) = error("Christoffel symbols not defined for this spacetime.")
set_christoffel!(Γ, position, ::AbstractSpacetime) = error("Christoffel symbols not defined for this spacetime.")

#Optional
set_metric_inverse!(::AbstractSpacetime) = error("Metric inverse not defined for this spacetime.")
event_horizon_radius(::AbstractSpacetime) = error("Event horizon radius not defined for this spacetime.")
circular_geodesic_angular_speed(position, ::AbstractSpacetime) = error("Circular geodesic angular speed not defined for this spacetime.")
allocate_christoffel_cache(::AbstractSpacetime) = nothing
set_christoffel!(Γ, position, spacetime::AbstractSpacetime, ::Nothing) = set_christoffel!(Γ, position, spacetime)

include("coordinatealias.jl")
include("general.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("numerical.jl")
include("RAR.jl")

"""
Computes the inverse of the given metric at the given position using a fast inversion
for 4x4 symmetric matrices.

Parameters:
- ginv: mutable array of size (4,4) to store the resulting inverse metric.
- position: tuple of four numbers representing a point in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) to store the metric evaluated at the given position.

Returns: nothing.
"""
function set_metric_inverse!(ginv, position, spacetime::AbstractSpacetime, g)
    set_metric!(g, position, spacetime)
    inverse_4x4_symmetric!(ginv, g)
    return nothing
end
