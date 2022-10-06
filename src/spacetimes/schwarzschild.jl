export SchwarzschildSpacetimeKerrSchildCoordinates 
export SchwarzschildSpacetimeSphericalCoordinates

@with_kw struct SchwarzschildSpacetimeKerrSchildCoordinates <: AnalyticSpacetime 

    M::Float64
    @assert M >= 0.0

end

@with_kw struct SchwarzschildSpacetimeSphericalCoordinates <: AnalyticSpacetime 

    M::Float64
    @assert M >= 0.0

end