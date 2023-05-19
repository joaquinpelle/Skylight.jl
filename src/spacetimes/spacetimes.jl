#Required
coordinates_topology(::AbstractSpacetime) = error("Coordinates topology not defined for this spacetime.")
set_metric!(::AbstractSpacetime) = error("Metric not defined for this spacetime.")
set_christoffel!(Γ, position, ::AbstractSpacetime, cache) = error("Christoffel symbols not defined for this spacetime.")
set_christoffel!(Γ, position, ::AbstractSpacetime) = error("Christoffel symbols not defined for this spacetime.")

#Optional
set_metric_inverse!(::AbstractSpacetime) = error("Metric inverse not defined for this spacetime.")
event_horizon_radius(::AbstractSpacetime) = error("Event horizon radius not defined for this spacetime.")
circular_geodesic_angular_speed(::AbstractSpacetime) = error("Circular geodesic angular speed not defined for this spacetime.")
allocate_christoffel_cache(::AbstractSpacetime) = nothing

include("coordinate_alias.jl")
include("general.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("numerical.jl")
include("RAR.jl")

set_christoffel!(Γ, position, spacetime::AbstractSpacetime, ::Nothing) = set_christoffel!(Γ, position, spacetime::AbstractSpacetime)