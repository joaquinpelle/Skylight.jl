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

struct JohannsenChristoffelCache <: ChristoffelCache end

allocate_christoffel_cache(spacetime::JohannsenSpacetimeBoyerLindquistCoordinates) = JohannsenChristoffelCache()

function set_christoffel!(Γ, point, spacetime::JohannsenSpacetimeBoyerLindquistCoordinates, cache) 
    
    t, r, θ, φ = point

    M = spacetime.M
    a = spacetime.a
    α13 = spacetime.α13
    α22 = spacetime.α22
    α52 = spacetime.α52
    ϵ3 = spacetime.ϵ3
    
    M2 = M^2
    M3 = M^3
    a2 = a^2
    
    r2 = r^2
    r3 = r^3
    r4 = r^4
    
    sinθ2 = sin(θ)^2
    sin2θ = sin(2*θ)
    
    Δ = r2 - 2*M*r + a2
    Σ = r2 + a2*cos(θ)^2 + ϵ3*M3/r
    
    A1 = 1 + α13*M3/r3
    A2 = 1 + α22*M2/r2
    A5 = 1 + α52*M2/r2
    
    Cp = (r2+a2)*A1-a2*A2*sinθ2
    C  = Cp^2 
    
    F = Δ-a2*A2^2*sinθ2
    G = (r2+a2)*A1*A2-Δ
    H = (r2+a2)^2*A1^2-a2*Δ*sinθ2
    
    
    gtt = -Σ*F/C
    gtφ = -a*G*Σ*sinθ2/C 

    grr = Σ/(Δ*A5)
    
    gθθ = Σ
    
    gφφ = Σ*sinθ2*H/C
    
    #The partial derivatives of the auxiliary functions
    
    ∂rΔ = 2*(r-M)
    ∂rΣ = 2*r-ϵ3*M3/r2
    ∂θΣ = -a2*sin2θ
    
    ∂rA1 = -3*α13*M3/r4
    ∂rA2 = -2*α22*M2/r3
    ∂rA5 = -2*α52*M2/r3
    
    
    ∂rCp = 2*r*A1+(r2+a2)*∂rA1-a2*∂rA2*sinθ2
    ∂θCp = -a2*A2*sin2θ
    
    ∂rF = ∂rΔ - 2*a2*sinθ2*A2*∂rA2
    ∂θF = -a2*A2^2*sin2θ
    
    ∂rG = 2*r*A1*A2+(r2+a2)*(∂rA1*A2+A1*∂rA2)-∂rΔ
    
    ∂rH = 4*(r2+a2)*r*A1^2+2*(r2+a2)^2*A1*∂rA1-a2*∂rΔ*sinθ2
    ∂θH = -a2*Δ*sin2θ
    
    #The partial derivatives of the metric
    
    ∂rgtt = -(Cp*(∂rΣ*F + Σ*∂rF) - 2*Σ*F*∂rCp)/(C*Cp)
    ∂θgtt = -(Cp*(∂θΣ*F + Σ*∂θF) - 2*Σ*F*∂θCp)/(C*Cp)
    
    ∂rgtφ = -a*sinθ2*(Cp*(∂rG*Σ + G*∂rΣ) - 2*G*Σ*∂rCp)/(C*Cp)
    ∂θgtφ = -a*G*(Cp*(∂θΣ*sinθ2 + Σ*sin2θ) - 2*Σ*sinθ2*∂θCp)/(C*Cp)
    
    ∂rgrr =(∂rΣ*Δ*A5 - Σ*(∂rΔ*A5+Δ*∂rA5))/(Δ*A5)^2
    ∂θgrr = ∂θΣ/(Δ*A5)
    
    ∂rgθθ = ∂rΣ
    ∂θgθθ = ∂θΣ
    
    ∂rgφφ = sinθ2*(Cp*(∂rΣ*H + Σ*∂rH) - 2*Σ*H*∂rCp)/(C*Cp)
    ∂θgφφ = (Cp*(sin2θ*Σ*H + ∂θΣ*sinθ2*H + ∂θH*Σ*sinθ2)-2*H*Σ*sinθ2*∂θCp)/(C*Cp)
    
    # Now we build the Christoffel symbols
    
    det = gtt*gφφ-gtφ^2
    
    Γ[1,1,2] = 0.5*(gφφ*∂rgtt-gtφ*∂rgtφ)/det
    Γ[1,2,1] = Γ[1,1,2]
    Γ[1,1,3] = 0.5*(gφφ*∂θgtt-gtφ*∂θgtφ)/det
    Γ[1,3,1] = Γ[1,1,3]
    Γ[1,2,4] = 0.5*(gφφ*∂rgtφ-gtφ*∂rgφφ)/det
    Γ[1,4,2] = Γ[1,2,4]
    Γ[1,3,4] = 0.5*(gφφ*∂θgtφ-gtφ*∂θgφφ)/det
    Γ[1,4,3] = Γ[1,3,4]
    
    Γ[2,1,1] = -0.5*∂rgtt/grr
    Γ[2,1,4] = -0.5*∂rgtφ/grr
    Γ[2,4,1] = Γ[2,1,4]
    Γ[2,2,2] = 0.5*∂rgrr/grr
    Γ[2,2,3] = 0.5*∂θgrr/grr
    Γ[2,3,2] = Γ[2,2,3]
    Γ[2,3,3] = -0.5*∂rgθθ/grr
    Γ[2,4,4] = -0.5*∂rgφφ/grr
    
    
    Γ[3,1,1] = -0.5*∂θgtt/gθθ
    Γ[3,1,4] = -0.5*∂θgtφ/gθθ
    Γ[3,4,1] = Γ[3,1,4]
    Γ[3,2,2] = -0.5*∂θgrr/gθθ
    Γ[3,2,3] = 0.5*∂rgθθ/gθθ
    Γ[3,3,2] = Γ[3,2,3]
    Γ[3,3,3] = 0.5*∂θgθθ/gθθ
    Γ[3,4,4] = -0.5*∂θgφφ/gθθ
    
    
    Γ[4,1,2] = -0.5*(gtφ*∂rgtt-gtt*∂rgtφ)/det
    Γ[4,2,1] = Γ[4,1,2]
    Γ[4,1,3] = -0.5*(gtφ*∂θgtt-gtt*∂θgtφ)/det
    Γ[4,3,1] = Γ[4,1,3]
    Γ[4,2,4] = -0.5*(gtφ*∂rgtφ-gtt*∂rgφφ)/det
    Γ[4,4,2] = Γ[4,2,4]
    Γ[4,3,4] = -0.5*(gtφ*∂θgtφ-gtt*∂θgφφ)/det
    Γ[4,4,3] = Γ[4,3,4]

    return nothing
    
end