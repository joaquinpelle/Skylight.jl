export ChargedWormholeSpacetimeSphericalCoordinates
export ChargedWormholeSpacetimeRegularCoordinates

@with_kw struct ChargedWormholeSpacetimeSphericalCoordinates <: AnalyticSpacetime
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end

@with_kw struct ChargedWormholeSpacetimeRegularCoordinates <: AnalyticSpacetime
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end
