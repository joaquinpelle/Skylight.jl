export KerrSpacetimeKerrSchildCoordinates
export KerrSpacetimeBoyerLindquistCoordinates

# Kerr Schild coordinates

@with_kw struct KerrSpacetimeKerrSchildCoordinates <: BlackHoleSpacetime

    M::Float64
    a::Float64

    @assert M >= 0.0
    @assert abs(a) <= M 

    #Metric cache
    l::Vector{Float64} = zeros(4)

end

@with_kw struct KerrKSChristoffelCache <: ChristoffelCache
    
    l::Vector{Float64} = zeros(4)
    dH::Array{Float64, 1} = zeros(4)
    dl::Array{Float64, 2} = zeros(4,4)
    D::Array{Float64, 3} = zeros(4,4,4)

end

coordinate_system_class(spacetime::KerrSpacetimeKerrSchildCoordinates) = CartesianClass()

event_horizon_radius(spacetime::KerrSpacetimeKerrSchildCoordinates) = spacetime.M*(1+sqrt(1-spacetime.a^2))

allocate_christoffel_cache(spacetime::KerrSpacetimeKerrSchildCoordinates) = KerrKSChristoffelCache()

function get_kerr_radius(position, spacetime::KerrSpacetimeKerrSchildCoordinates)
    
    t, x, y, z = position
    ρ2 = x^2 + y^2 + z^2
    a2 = spacetime.a^2
    r2 = 0.5 * (ρ2 - a2) + sqrt(0.25 * (ρ2 - a2)^2 + a2 * z^2)

    return sqrt(r2) 
end

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
    
    return nothing
    
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
    
    return nothing
    
end

function set_christoffel!(Γ, position, spacetime::KerrSpacetimeKerrSchildCoordinates, cache::KerrKSChristoffelCache)

    t, x, y, z = position
    M = spacetime.M
    a = spacetime.a

    a2 = a^2
    rho2 = x^2 .+ y^2 .+ z^2
    r2 = 0.5*(rho2-a2) .+ sqrt(0.25*(rho2-a2)^2+a2*z^2)
    r  = sqrt(r2)
    r3 = r2*r
    r4 = r2*r2

    #Derivatives of r(x,y,z)
    dr_dx = x*r3*(r2+a2)/(a2*z^2*(2*r2+a2)+r4*rho2)
    dr_dy = y*r3*(r2+a2)/(a2*z^2*(2*r2+a2)+r4*rho2)
    dr_dz = z*r*(r2+a2)^2/(a2*z^2*(2*r2+a2)+r4*rho2)

    #The scalar function and the null vector of the metric

    H = M*r3/(r4+a2*z^2)

    l = cache.l
    l[1] = 1.
    l[2] = (r*x + a*y)/(r2 + a2)
    l[3] = (r*y - a*x)/(r2 + a2)
    l[4] = z/r

    #Derivatives of the null vectors ( dl[a,b]=d_a l[b]). Indexes are down.

    dl = cache.dl
    fill!(dl, 0)

    dl[2,2] =(dr_dx*(x-2 * r*l[2])+r)/(r2+a2)
    dl[2,3] =(dr_dx*(y-2 * r*l[3])-a)/(r2+a2)
    dl[2,4] =-z/r2*dr_dx

    dl[3,2] =(dr_dy*(x-2 *r*l[2])+a)/(r2+a2)
    dl[3,3] =(dr_dy*(y-2 *r*l[3])+r)/(r2+a2)
    dl[3,4] =-z/r2*dr_dy

    dl[4,2] =dr_dz*(x-2 *r*l[2])/(r2+a2)
    dl[4,3] =dr_dz*(y-2 *r*l[3])/(r2+a2)
    dl[4,4] =1.0/r-z/r2*dr_dz


    #Derivatives of the scalar function H (dH[a]=d_a H). Index is down.

    dH = cache.dH
    fill!(dH, 0)

    dH[2] = -M*r2*(r4-3*a2*z^2)/(r4+a2*z^2)^2*dr_dx
    dH[3] = -M*r2*(r4-3*a2*z^2)/(r4+a2*z^2)^2*dr_dy
    dH[4] = -M*r2*(2*a2*r*z+(r4-3*a2*z^2)*dr_dz)/(r4+a2*z^2)^2

    # Directional derivative of H in the direction of the null vector l  (l^a d_a H)
    l_dH = -M*r2*(r4-a2*z^2)/(r4+a2*z^2)^2

    # Tensor product of the null vector derivatives with the null vector. dlablc[a,b,c]= dl[a,b]*l[c] 
    # Derivatives of the products H*la*lb:  D[a,b,c]= d_a (H*lb*lc) (The order of fors gives the order of indexes)
    # This computation is equivalent to D[a,b,c]=dH[a]*l[b]*l[c]+H*dl[a,b]*l[c]+H*dl[a,c]*l[b]

    D = cache.D
    for i = 1:4
        for j = 1:4
            for k = 1:4
                D[i,j,k] = dH[i] * l[j] * l[k] + H * dl[i,j] * l[k] + H * dl[i,k] * l[j]
            end
        end
    end

    #Christoffel symbols

    for i = 1:4
        sign = i == 1 ? -1 : 1
        for j = 1:4
            for k = 1:4
                Γ[i,j,k] = sign * (D[j,k,i] + D[k,j,i] - D[i,j,k] + 2 * H * l_dH * l[i]*l[j]*l[k])
            end
        end
    end

    return nothing

end

# Boyer Lindquist coordinates

@with_kw struct KerrSpacetimeBoyerLindquistCoordinates <: BlackHoleSpacetime 

    M::Float64
    a::Float64

    @assert M >= 0.0
    @assert abs(a) <= M 

end   


coordinate_system_class(spacetime::KerrSpacetimeBoyerLindquistCoordinates) = SphericalClass()

event_horizon_radius(spacetime::KerrSpacetimeBoyerLindquistCoordinates) = spacetime.M*(1+sqrt(1-spacetime.a^2))

function set_metric!(g, point, spacetime::KerrSpacetimeBoyerLindquistCoordinates)
    
    t,r,θ,φ = point

    M = spacetime.M
    a = spacetime.a   
    
    a2 = a^2
    r2 = r^2
    senθ2 = sin(θ)^2
    
    Σ = r2+a2*cos(θ)^2
    Δ = r2-2*M*r+a2
     
    fill!(g,0.0)
    
    g[1,1] = -1.0 + 2*M*r/Σ
    g[1,4] = -2*M*a*r*senθ2/Σ
    
    g[2,2] = Σ/Δ 
    
    g[3,3] = Σ
    
    g[4,1] = g[1,4]
    g[4,4] = senθ2*(r2+a2+2*M*a2*r*senθ2/Σ)
    
    return nothing

end

function set_metric_inverse!(ginv, point, spacetime::KerrSpacetimeBoyerLindquistCoordinates)
    
    t,r,θ,φ = point

    M = spacetime.M
    a = spacetime.a   
    
    a2 = a^2
    r2 = r^2
    senθ2 = sin(θ)^2
    
    Σ = r2+a2*cos(θ)^2
    Δ = r2-2*M*r+a2
     
    fill!(ginv,0.0)
    
    ctt = senθ2*(r2+a2+2*M*a2*r*senθ2/Σ)
    cφφ = -1.0 + 2*M*r/Σ
    ctφ = 2*M*a*r*senθ2/Σ

    det = ctt*cφφ - ctφ^2

    ginv[1,1] = ctt/det
    ginv[1,4] = ctφ/det
    
    ginv[2,2] = Δ/Σ 
    
    ginv[3,3] = 1/Σ
    
    ginv[4,1] = ginv[1,4]
    ginv[4,4] = cφφ/det
    
    return nothing

end