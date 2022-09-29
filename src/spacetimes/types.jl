abstract type Spacetime end

abstract type AnalyticSpacetime <: Spacetime end
abstract type NumericalSpacetime <: Spacetime end

struct FullNumericalSpacetime <: NumericalSpacetime end
struct StaticSphericallySymmetricSpacetime <: NumericalSpacetime end


abstract type AsymptoticCoordinateSystem end

struct CartesianCoordinates <: AsymptoticCoordinateSystem end
struct SphericalCoordinates <: AsymptoticCoordinateSystem end


abstract type SpacetimeParameters end

