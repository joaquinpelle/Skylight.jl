export ChargedWormholeSpacetimeSphericalCoordinates
export ChargedWormholeSpacetimeRegularCoordinates

@with_kw struct ChargedWormholeSpacetimeSphericalCoordinates <: WormholeSpacetime
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end

coordinate_system_class(spacetime::ChargedWormholeSpacetimeSphericalCoordinates) = SphericalClass()

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

@with_kw struct ChargedWormholeSpacetimeRegularCoordinates <: WormholeSpacetime
    
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0

end

function get_wormhole_radius(l, spacetime::WormholeSpacetime)

    b0 = spacetime.b0
    Q = spacetime.Q

    return sqrt(l^2+b0^2-Q^2)

end

coordinate_system_class(spacetime::ChargedWormholeSpacetimeRegularCoordinates) = SphericalClass()

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

struct ChargedWormholeChristoffelCache <: ChristoffelCache end

allocate_christoffel_cache(spacetime::ChargedWormholeSpacetimeRegularCoordinates) = ChargedWormholeChristoffelCache()

function set_christoffel!(Γ,position,spacetime::ChargedWormholeSpacetimeRegularCoordinates,cache::ChargedWormholeChristoffelCache)

    #Spacetime coordinates
    t, l, θ, φ = position

    b0 = spacetime.b0
    Q = spacetime.Q
    
    r = sqrt(l^2+b0^2-Q^2)

    gtt = -(1+Q^2/r^2)
    gθθ = r^2
    gφφ = r^2*sin(θ)^2

    dr_dl = l/r

    dlgtt = 2*Q^2/r^3*dr_dl
    dlgθθ = 2*l
    dlgφφ = 2*l*sin(θ)^2
    dθgφφ = r^2*2*sin(θ)*cos(θ)

    #The Christoffel symbols

    Γ[1,1,2] = 0.5*dlgtt/gtt
    Γ[1,2,1] = Γ[1,1,2]

    Γ[2,1,1] = -dlgtt/2
    Γ[2,3,3] = -dlgθθ/2
    Γ[2,4,4] = -dlgφφ/2

    Γ[3,2,3] = 0.5*dlgθθ/gθθ
    Γ[3,4,4] = -sin(θ)*cos(θ)
    Γ[3,3,2] = Γ[3,2,3]

    Γ[4,2,4] = 0.5*dlgφφ/gφφ
    Γ[4,3,4] = cos(θ)/sin(θ)
    Γ[4,4,2] = Γ[4,2,4]
    Γ[4,4,3] = Γ[4,3,4]

    return nothing

end
