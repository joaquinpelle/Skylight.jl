abstract type Spacetime end

abstract type AnalyticSpacetime <: Spacetime end
abstract type NumericalSpacetime <: Spacetime end

abstract type FlatSpacetime <: AnalyticSpacetime end
abstract type BlackHoleSpacetime <: AnalyticSpacetime end
abstract type WormholeSpacetime <: AnalyticSpacetime end

abstract type CoordinateSystemClass end

struct CartesianClass <: CoordinateSystemClass end
struct SphericalClass <: CoordinateSystemClass end

abstract type ChristoffelCache end

#Required
coordinate_system_class(spacetime::Spacetime) = error("Coordinate system not defined for this spacetime.")
set_metric!(spacetime::Spacetime) = error("Metric not defined for this spacetime.")
allocate_christoffel_cache(spacetime::Spacetime) = error("Christoffel cache not defined for this spacetime.")
set_christoffel!(Γ, position, spacetime::Spacetime, cache) = error("Christoffel symbols not defined for this spacetime.")

#Optional
set_metric_inverse!(spacetime::Spacetime) = error("Metric inverse not defined for this spacetime.")
event_horizon_radius(spacetime::Spacetime) = error("Event horizon radius not defined for this spacetime.")
circular_geodesic_angular_speed(spacetime::Spacetime) = error("Circular geodesic angular speed not defined for this spacetime.")


include("coordinate_alias.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("fullnumerical.jl")
include("RAR.jl")