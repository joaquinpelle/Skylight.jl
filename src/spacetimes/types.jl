abstract type AbstractSpacetime end
abstract type AbstractAutoDiffSpacetime <: AbstractSpacetime end
abstract type AbstractChristoffelCache end
abstract type AbstractCoordinatesTopology end

struct CartesianTopology <: AbstractCoordinatesTopology end
struct SphericalTopology <: AbstractCoordinatesTopology end

abstract type AbstractRotationSense end

struct ProgradeRotation <: AbstractRotationSense end
struct RetrogradeRotation <: AbstractRotationSense end

abstract type AbstractSpacetimeSymmetry end

abstract type Stationarity <: AbstractSpacetimeSymmetry end

struct IsStationary <: Stationarity end 
struct IsNotStationary <: Stationarity end 

abstract type SphericalSymmetry <: AbstractSpacetimeSymmetry end

struct IsSphericallySymmetric <: SphericalSymmetry end
struct IsNotSphericallySymmetric <: SphericalSymmetry end

abstract type AxialSymmetry <: AbstractSpacetimeSymmetry end

struct IsAxiallySymmetric <: AxialSymmetry end
struct IsNotAxiallySymmetric <: AxialSymmetry end

