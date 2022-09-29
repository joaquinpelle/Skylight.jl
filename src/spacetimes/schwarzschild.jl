export SchwarzschildSpacetimeParameters
export SchwarzschildSpacetimeKerrSchildCoordinates, SchwarzschildSpacetimeSphericalCoordinates

@with_kw struct SchwarzschildSpacetimeParameters <: SpacetimeParameters
    
    M::Float64
    @assert M >= 0.0

end

# Kerr-Schild coordinates

struct SchwarzschildSpacetimeKerrSchildCoordinates <: AnalyticSpacetime end

# Spherical coordinates

struct SchwarzschildSpacetimeSphericalCoordinates <: AnalyticSpacetime end