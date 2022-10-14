export ChargedWormholeSpacetimeSphericalCoordinates
export ChargedWormholeSpacetimeRegularCoordinates

@with_kw struct ChargedWormholeSpacetimeSphericalCoordinates <: AnalyticSpacetime
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end

coordinate_system_kind(spacetime::ChargedWormholeSpacetimeSphericalCoordinates) = SphericalKind()

function set_metric!(g, point, spacetime::ChargedWormholeSpacetimeSphericalCoordinates)
        
    t, r, θ, φ = point

    b0 = spacetime.b0
    Q = spacetime.Q

    b = b0^2/r

    gtt = -(1+Q^2/r^2)
    grr = 1/(1-b/r+Q^2/r^2)    
    
    g[1,1]= gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= grr 
    g[2,3]= 0.
    g[2,4]= 0.
    g[3,1]= 0.
    g[3,2]= 0.
    g[3,3]= r^2
    g[3,4]= 0.
    g[4,1]= 0.
    g[4,2]= 0.
    g[4,3]= 0.
    g[4,4]= r^2*sin(θ)^2
    
    return nothing

end

function set_metric_inverse!(g, point, spacetime::ChargedWormholeSpacetimeSphericalCoordinates)
        
    t, r, θ, φ = point

    b0 = spacetime.b0
    Q = spacetime.Q

    b = b0^2/r

    gtt = -(1+Q^2/r^2)
    grr = 1/(1-b/r+Q^2/r^2)

    g[1,1]= 1.0/gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= 1.0/grr 
    g[2,3]= 0.
    g[2,4]= 0.
    g[3,1]= 0.
    g[3,2]= 0.
    g[3,3]= 1.0/r^2
    g[3,4]= 0.
    g[4,1]= 0.
    g[4,2]= 0.
    g[4,3]= 0.
    g[4,4]= 1.0/(r^2*sin(θ)^2)
    
    return nothing

end

@with_kw struct ChargedWormholeSpacetimeRegularCoordinates <: AnalyticSpacetime
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end

coordinate_system_kind(spacetime::ChargedWormholeSpacetimeRegularCoordinates) = SphericalKind()

function set_metric!(g, point, spacetime::ChargedWormholeSpacetimeRegularCoordinates)
        
    t, l, θ, φ = point

    b0 = spacetime.b0
    Q = spacetime.Q

    r = sqrt(l^2+b0^2-Q^2)

    gtt = -(1+Q^2/r^2)
    gll = 1.0
    
    g[1,1]= gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= gll 
    g[2,3]= 0.
    g[2,4]= 0.
    g[3,1]= 0.
    g[3,2]= 0.
    g[3,3]= r^2
    g[3,4]= 0.
    g[4,1]= 0.
    g[4,2]= 0.
    g[4,3]= 0.
    g[4,4]= r^2*sin(θ)^2
    
    return nothing

end

function set_metric_inverse!(g, point, spacetime::ChargedWormholeSpacetimeRegularCoordinates)
        
    t, l, θ, φ = point

    b0 = spacetime.b0
    Q = spacetime.Q

    r = sqrt(l^2+b0^2-Q^2)

    gtt = -(1+Q^2/r^2)
    gll = 1.0
    
    g[1,1]= 1.0/gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= 1.0/gll 
    g[2,3]= 0.
    g[2,4]= 0.
    g[3,1]= 0.
    g[3,2]= 0.
    g[3,3]= 1.0/r^2
    g[3,4]= 0.
    g[4,1]= 0.
    g[4,2]= 0.
    g[4,3]= 0.
    g[4,4]= 1.0/(r^2*sin(θ)^2)
    
    return nothing
    
end