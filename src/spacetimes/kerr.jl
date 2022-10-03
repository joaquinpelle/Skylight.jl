export KerrSpacetimeParameters 
export KerrSpacetimeKerrSchildCoordinates, KerrSpacetimeBoyerLindquistCoordinates

@with_kw struct KerrSpacetimeParameters <: SpacetimeParameters
    
    M::Float64
    a::Float64

    @assert M >= 0.0
    @assert abs(a) <= M

end

# Kerr Schild coordinates

@with_kw struct KerrSpacetimeKerrSchildCoordinates{T<:Function} <: AnalyticSpacetime

    parameters::KerrSpacetimeParameters 
    coordinate_system_kind::CartesianKind = CartesianKind()
    metric!::T = kerr_metric_kerr_schild_coordinates!

end

function kerr_metric_kerr_schild_coordinates!(g, q, par::KerrSpacetimeParameters)

    """ 
    g: container for the metric 
    q: spacetime position
    """

    M = par.M
    a = par.a
    
    t, x, y, z = q
    ρ2 = x^2 + y^2 + z^2
    a2 = a^2
    r2 = 0.5 * (ρ2 - a2) + sqrt(0.25 * (ρ2 - a2)^2 + a2 * z^2)
    r = sqrt(r2)
    H2 = 2. * M * r / (r2 + a2 * z^2 / r2)
    
    l = zeros(4)
    l[1] = 1.
    l[2] = (r*x + a*y)/(r2 + a2)
    l[3] = (r*y - a*x)/(r2 + a2)
    l[4] = z/r
    
    g[1,1]=-1. + H2 * l[1]*l[1]
    g[1,2]= 0. + H2 * l[1]*l[2]
    g[1,3]= 0. + H2 * l[1]*l[3]
    g[1,4]= 0. + H2 * l[1]*l[4]
    g[2,1]= g[1,2]
    g[2,2]= 1. + H2 * l[2]*l[2]
    g[2,3]= 0. + H2 * l[2]*l[3]
    g[2,4]= 0. + H2 * l[2]*l[4]
    g[3,1]= g[1,3]
    g[3,2]= g[2,3]
    g[3,3]= 1. + H2 * l[3]*l[3]
    g[3,4]= 0. + H2 * l[3]*l[4]
    g[4,1]= g[1,4]
    g[4,2]= g[2,4]
    g[4,3]= g[3,4]
    g[4,4]= 1. + H2 * l[4]*l[4]
    
    return g
    
end

# Boyer Lindquist coordinates

struct KerrSpacetimeBoyerLindquistCoordinates <: AnalyticSpacetime end   
