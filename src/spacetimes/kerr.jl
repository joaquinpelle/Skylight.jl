export KerrSpacetimeKerrSchildCoordinates
export KerrSpacetimeBoyerLindquistCoordinates

# Kerr Schild coordinates

@with_kw struct KerrSpacetimeKerrSchildCoordinates <: AnalyticSpacetime

    M::Float64
    a::Float64

    @assert M >= 0.0
    @assert abs(a) <= M 

    #Cache
    l::Vector{Float64} = zeros(4)

end

coordinate_system_kind(spacetime::KerrSpacetimeKerrSchildCoordinates) = CartesianKind()

function set_metric!(g, position, spacetime::KerrSpacetimeKerrSchildCoordinates)

    """ 
    g: container for the metric 
    q: spacetime position
    """

    M = spacetime.M
    a = spacetime.a
    
    t, x, y, z = position
    ρ2 = x^2 + y^2 + z^2
    a2 = a^2
    r2 = 0.5 * (ρ2 - a2) + sqrt(0.25 * (ρ2 - a2)^2 + a2 * z^2)
    r = sqrt(r2)
    H2 = 2. * M * r / (r2 + a2 * z^2 / r2)
    
    l = spacetime.l
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

function set_metric_inverse!(g, position, spacetime::KerrSpacetimeKerrSchildCoordinates)

    M = spacetime.M
    a = spacetime.a
    
    t, x, y, z = position
    ρ2 = x^2 + y^2 + z^2
    a2 = a^2
    r2 = 0.5 * (ρ2 - a2) + sqrt(0.25 * (ρ2 - a2)^2 + a2 * z^2)
    r = sqrt(r2)
    H2 = 2. * M * r / (r2 + a2 * z^2 / r2)
    
    l = spacetime.l
    l[1] = -1.0
    l[2] = (r*x + a*y)/(r2 + a2)
    l[3] = (r*y - a*x)/(r2 + a2)
    l[4] = z/r
    
    g[1,1]=-1. - H2 * l[1]*l[1]
    g[1,2]= 0. - H2 * l[1]*l[2]
    g[1,3]= 0. - H2 * l[1]*l[3]
    g[1,4]= 0. - H2 * l[1]*l[4]
    g[2,1]= g[1,2]
    g[2,2]= 1. - H2 * l[2]*l[2]
    g[2,3]= 0. - H2 * l[2]*l[3]
    g[2,4]= 0. - H2 * l[2]*l[4]
    g[3,1]= g[1,3]
    g[3,2]= g[2,3]
    g[3,3]= 1. - H2 * l[3]*l[3]
    g[3,4]= 0. - H2 * l[3]*l[4]
    g[4,1]= g[1,4]
    g[4,2]= g[2,4]
    g[4,3]= g[3,4]
    g[4,4]= 1. - H2 * l[4]*l[4]
    
    return g
    
end



# Boyer Lindquist coordinates

struct KerrSpacetimeBoyerLindquistCoordinates <: AnalyticSpacetime end   
