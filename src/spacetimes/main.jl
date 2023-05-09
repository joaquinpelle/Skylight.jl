abstract type AbstractSpacetime end

abstract type AnalyticSpacetime <: AbstractSpacetime end
abstract type NumericalSpacetime <: AbstractSpacetime end

abstract type FlatSpacetime <: AnalyticSpacetime end
abstract type BlackHoleSpacetime <: AnalyticSpacetime end
abstract type WormholeSpacetime <: AnalyticSpacetime end

abstract type CoordinateSystemClass end

struct CartesianClass <: CoordinateSystemClass end
struct SphericalClass <: CoordinateSystemClass end

abstract type ChristoffelCache end

#Required
coordinate_system_class(spacetime::AbstractSpacetime) = error("Coordinate system not defined for this spacetime.")
set_metric!(spacetime::AbstractSpacetime) = error("Metric not defined for this spacetime.")
allocate_christoffel_cache(spacetime::AbstractSpacetime) = error("Christoffel cache not defined for this spacetime.")
set_christoffel!(Γ, position, spacetime::AbstractSpacetime, cache) = error("Christoffel symbols not defined for this spacetime.")

#Optional
set_metric_inverse!(spacetime::AbstractSpacetime) = error("Metric inverse not defined for this spacetime.")
event_horizon_radius(spacetime::AbstractSpacetime) = error("Event horizon radius not defined for this spacetime.")
circular_geodesic_angular_speed(spacetime::AbstractSpacetime) = error("Circular geodesic angular speed not defined for this spacetime.")


include("coordinate_alias.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("fullnumerical.jl")
include("RAR.jl")

export CartesianClass, SphericalClass
export coordinate_system_class, set_metric!, allocate_christoffel_cache, set_christoffel!