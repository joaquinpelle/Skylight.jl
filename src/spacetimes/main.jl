abstract type Spacetime end

abstract type AnalyticSpacetime <: Spacetime end
abstract type NumericalSpacetime <: Spacetime end

abstract type CoordinateSystemKind end

struct CartesianKind <: CoordinateSystemKind end
struct SphericalKind <: CoordinateSystemKind end

abstract type ChristoffelCache end

include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("fullnumerical.jl")
include("staticsphericallysymmetric.jl")