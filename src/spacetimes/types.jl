abstract type AbstractSpacetime end
abstract type AbstractAutoDiffSpacetime <: AbstractSpacetime end
abstract type AbstractChristoffelCache end
abstract type AbstractCoordinatesTopology end

struct CartesianTopology <: AbstractCoordinatesTopology end
struct SphericalTopology <: AbstractCoordinatesTopology end

@with_kw mutable struct AutoDiffChristoffelCache <: AbstractChristoffelCache
    g::Array{Float64,2} = zeros(4,4)
    ginv::Array{Float64,2} = zeros(4,4)
    ∂g::Array{Float64,3} = zeros(4,4,4)
end
