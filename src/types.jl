abstract type Spacetime end

abstract type AnalyticSpacetime <: Spacetime end
abstract type NumericalSpacetime <: Spacetime end

struct FullNumericalSpacetime <: NumericalSpacetime end
struct StaticSphericallySymmetricSpacetime <: NumericalSpacetime end


abstract type SpacetimeParameters end


abstract type CoordinateSystemKind end

struct CartesianKind <: CoordinateSystemKind end
struct SphericalKind <: CoordinateSystemKind end



