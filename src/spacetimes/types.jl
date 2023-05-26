abstract type AbstractSpacetime end
abstract type AbstractAutoDiffSpacetime <: AbstractSpacetime end
abstract type AbstractChristoffelCache end
abstract type AbstractCoordinatesTopology end

struct CartesianTopology <: AbstractCoordinatesTopology end
struct SphericalTopology <: AbstractCoordinatesTopology end

abstract type AbstractRotationSense end

struct ProgradeRotation <: AbstractRotationSense end
struct RetrogradeRotation <: AbstractRotationSense end
