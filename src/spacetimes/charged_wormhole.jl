export ChargedWormholeSpacetimeParameters
export ChargedWormholeSpacetimeSphericalCoordinates, ChargedWormholeSpacetimeSphericalCoordinates

@with_kw struct ChargedWormholeSpacetimeParameters <: SpacetimeParameters
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end

struct ChargedWormholeSpacetimeSphericalCoordinates <: AnalyticSpacetime end
struct ChargedWormholeSpacetimeRegularCoordinates <: AnalyticSpacetime end