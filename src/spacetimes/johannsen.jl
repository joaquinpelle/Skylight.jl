export JohannsenSpacetimeBoyerLindquistCoordinates

@with_kw struct JohannsenSpacetimeBoyerLindquistCoordinates <: BlackHoleSpacetime
    
    """Johannsen (2013) spacetime to lowest order in the deformation parameters from the Kerr metric

    Reference: https://arxiv.org/pdf/1501.02809.pdf

    q: spacetime position in the background coordinates
    M: the mass of the spacetime
    a: the black hole spin
    α13, α22, α52, ϵ3: the lowest order deformation parameters

    """

    M::Float64
    a::Float64
    α13::Float64
    α22::Float64
    α52::Float64
    ϵ3::Float64

    @assert M >= 0.0
    @assert abs(a) <= M    

    factor = (M+sqrt(M^2-a^2))/M
    
    @assert  α13 > -factor^3
    @assert  α22 > -factor^2
    @assert  α52 > -factor^2
    @assert  ϵ3 > -factor^3
        
end

event_horizon_radius(spacetime::JohannsenSpacetimeBoyerLindquistCoordinates) = spacetime.M*(1+sqrt(1-spacetime.a^2))

coordinate_system_class(spacetime::JohannsenSpacetimeBoyerLindquistCoordinates) = SphericalClass()

function set_metric!(g,point,spacetime::JohannsenSpacetimeBoyerLindquistCoordinates)

    t, r, θ, φ = point

    M = spacetime.M
    a = spacetime.a
    α13 = spacetime.α13
    α22 = spacetime.α22
    α52 = spacetime.α52
    ϵ3 = spacetime.ϵ3

    r2 = r^2
    a2 = a^2
    sinθ2 = sin(θ)^2

    Δ = r2 - 2*M*r + a2
    Σ = r2 + a2*cos(θ)^2 + ϵ3*M^3/r

    A1 = 1 + α13*(M/r)^3
    A2 = 1 + α22*(M/r)^2
    A5 = 1 + α52*(M/r)^2

    C = ((r2+a2)*A1-a2*A2*sinθ2)^2

    fill!(g,0.0)

    g[1,1] = -Σ*(Δ-a2*A2^2*sinθ2)/C
    g[1,4] = -a*((r2+a2)*A1*A2-Δ)*Σ*sinθ2/C 

    g[2,2] = Σ/(Δ*A5)

    g[3,3] = Σ

    g[4,1] = g[1,4]
    g[4,4] = Σ*sinθ2*((r2+a2)^2*A1^2-a2*Δ*sinθ2)/C

    return nothing

end

function set_metric_inverse!(g, point, spacetime::JohannsenSpacetimeBoyerLindquistCoordinates)

    t, r, θ, φ = point

    M = spacetime.M
    a = spacetime.a
    α13 = spacetime.α13
    α22 = spacetime.α22
    α52 = spacetime.α52
    ϵ3 = spacetime.ϵ3

    r2 = r^2
    a2 = a^2
    sinθ2 = sin(θ)^2

    Δ = r2 - 2*M*r + a2
    Σ = r2 + a2*cos(θ)^2 + ϵ3*M^3/r

    A1 = 1 + α13*(M/r)^3
    A2 = 1 + α22*(M/r)^2
    A5 = 1 + α52*(M/r)^2

    C = ((r2+a2)*A1-a2*A2*sinθ2)^2

    # Finally, the metric

    gtt = -Σ*(Δ-a2*A2^2*sinθ2)/C
    gtφ = -a*((r2+a2)*A1*A2-Δ)*Σ*sinθ2/C 
    gφφ = Σ*sinθ2*((r2+a2)^2*A1^2-a2*Δ*sinθ2)/C

    det = gtt*gφφ-gtφ^2

    fill!(g,0)

    g[1,1] = gφφ/det 
    g[1,4] = -gtφ/det 

    g[2,2] = (Δ*A5)/Σ

    g[3,3] = 1.0/Σ

    g[4,1] = g[1,4]
    g[4,4] = gtt/det 

    return nothing

end