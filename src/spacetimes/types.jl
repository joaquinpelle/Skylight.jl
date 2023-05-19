abstract type AbstractSpacetime end

abstract type AbstractCoordinateTopology end

struct CartesianTopology <: AbstractCoordinateTopology end
struct SphericalTopology <: AbstractCoordinateTopology end

abstract type ChristoffelCache end