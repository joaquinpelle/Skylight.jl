export JohannsenSpacetimeBoyerLindquistCoordinates

@with_kw struct JohannsenSpacetimeBoyerLindquistCoordinates <: AnalyticSpacetime
    
    M::Float64
    a::Float64
    α13::Float64
    α22::Float64
    α52::Float64
    ϵ3::Float64

    factor = (M+sqrt(M^2-a^2))/M
    @assert M >= 0.0
    @assert abs(a) <= M    
    @assert  α13 > -factor^3
    @assert  α22 > -factor^2
    @assert  α52 > -factor^2
    @assert  ϵ3 > -factor^3
        
end
