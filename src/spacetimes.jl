abstract type Spacetime end

abstract type AnalyticSpacetime <: Spacetime end
abstract type NumericalSpacetime <: Spacetime end

abstract type SpacetimeParameters end

include("spacetimes/minkowski.jl")
include("spacetimes/schwarzschild.jl")
include("spacetimes/kerr.jl")
include("spacetimes/johannsen.jl")
include("spacetimes/chargedwormhole.jl")
include("spacetimes/fullnumerical.jl")
include("spacetimes/staticsphericallysymmetric.jl")