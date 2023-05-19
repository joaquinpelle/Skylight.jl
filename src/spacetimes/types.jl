abstract type AbstractSpacetime end
abstract type AbstractChristoffelCache end
abstract type AbstractCoordinatesTopology end

struct CartesianTopology <: AbstractCoordinatesTopology end
struct SphericalTopology <: AbstractCoordinatesTopology end
