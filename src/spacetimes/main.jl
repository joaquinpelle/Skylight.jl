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

include("coordinate_alias.jl")
include("minkowski.jl")
include("schwarzschild.jl")
include("kerr.jl")
include("johannsen.jl")
include("chargedwormhole.jl")
include("bosonstar.jl")
include("fullnumerical.jl")
include("staticsphericallysymmetric.jl")