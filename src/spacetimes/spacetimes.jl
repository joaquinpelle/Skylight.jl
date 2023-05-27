#Required
coordinates_topology(spacetime) = error("Coordinates topology not defined for this spacetime.")
set_metric!(g, position, spacetime) = error("Metric not defined for this spacetime.")

#Optional
# set_metric_inverse!(ginv, position, ::AbstractSpacetime) = error("Metric inverse not defined for this spacetime.")
# volume_element(position, ::AbstractSpacetime, g) = error("Volume element not defined for this spacetime.")
radius(position, spacetime) = error("Radius not defined for this spacetime.")
event_horizon_radius(spacetime) = error("Event horizon radius not defined for this spacetime.")
circular_geodesic_angular_speed(position, spacetime, rotation_sense) = error("Circular geodesic angular speed not defined for this spacetime.")
#custom set_christoffel
set_christoffel!(Γ, position, spacetime::AbstractSpacetime, ::Nothing) = set_christoffel!(Γ, position, spacetime)

#By default we set non-stationarity and non-spherical symmetry  
is_stationary(::AbstractSpacetime) = IsNotStationary()
is_spherically_symmetric(::AbstractSpacetime) = IsNotSphericallySymmetric()
#For z symmetry we check spherical symmetry by default first. Thus if a spacetime is declared spherically symmetric it's automatically axially symmetric.
is_axially_symmetric(spacetime::AbstractSpacetime) = isa(is_spherically_symmetric(spacetime), IsSphericallySymmetric) ? IsAxiallySymmetric() : IsNotAxiallySymmetric()

include("coordinatealias.jl")
include("autodiff.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("RAR.jl")
include("numerical.jl")
include("constantsmotion.jl")

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
"""
Computes the volume element (square root of minus the determinant of the metric) at a given position 
using a fast determinant for 4x4 symmetric matrices.

Parameters:
- position: tuple of four numbers representing a point in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) to store the metric evaluated at the given position.

Returns: the volume element.
"""
function volume_element(position, spacetime::AbstractSpacetime, g)
    set_metric!(g, position, spacetime)
    return sqrt(-det_4x4_symmetric(g))
end

@inline sign(::ProgradeRotation) = 1.0
@inline sign(::RetrogradeRotation) = -1.0