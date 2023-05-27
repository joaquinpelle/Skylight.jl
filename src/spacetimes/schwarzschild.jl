abstract type AbstractSchwarzschildSpacetime <: AbstractSpacetime end

#Kerr-Schild coordinates

@with_kw struct SchwarzschildSpacetimeKerrSchildCoordinates <: AbstractSchwarzschildSpacetime 
    M::Float64
    @assert M >= 0.0

    #Cache
    l::Vector{Float64}=zeros(4)
end

coordinates_topology(::SchwarzschildSpacetimeKerrSchildCoordinates) = CartesianTopology()
radius(position,::SchwarzschildSpacetimeKerrSchildCoordinates) = sqrt(position[2]^2 + position[3]^2 + position[4]^2)

function set_metric!(g, position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates)
    
    M = spacetime.M
    
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    H = 2M/r
    
    l = spacetime.l

    l[1] = 1.0
    l[2] = x/r
    l[3] = y/r
    l[4] = z/r
    
    g[1,1]=-1. + H*l[1]*l[1]
    g[1,2]= 0. + H*l[1]*l[2]
    g[1,3]= 0. + H*l[1]*l[3]
    g[1,4]= 0. + H*l[1]*l[4]
    g[2,1]= g[1,2]
    g[2,2]= 1. + H*l[2]*l[2]
    g[2,3]= 0. + H*l[2]*l[3]
    g[2,4]= 0. + H*l[2]*l[4]
    g[3,1]= g[1,3]
    g[3,2]= g[2,3]
    g[3,3]= 1. + H*l[3]*l[3]
    g[3,4]= 0. + H*l[3]*l[4]
    g[4,1]= g[1,4]
    g[4,2]= g[2,4]
    g[4,3]= g[3,4]
    g[4,4]= 1. + H*l[4]*l[4]
    
    return nothing
end

function set_metric_inverse!(g, position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates)
    
    M = spacetime.M
    
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    H = 2M/r
    
    l = spacetime.l
    
    l[1] =-1.
    l[2] = x/r
    l[3] = y/r
    l[4] = z/r
    
    g[1,1]=-1. - H
    g[1,2]= 0. - H*l[1]*l[2]
    g[1,3]= 0. - H*l[1]*l[3]
    g[1,4]= 0. - H*l[1]*l[4]
    g[2,1]= g[1,2]
    g[2,2]= 1. - H*l[2]*l[2]
    g[2,3]= 0. - H*l[2]*l[3]
    g[2,4]= 0. - H*l[2]*l[4]
    g[3,1]= g[1,3]
    g[3,2]= g[2,3]
    g[3,3]= 1. - H*l[3]*l[3]
    g[3,4]= 0. - H*l[3]*l[4]
    g[4,1]= g[1,4]
    g[4,2]= g[2,4]
    g[4,3]= g[3,4]
    g[4,4]= 1. - H*l[4]*l[4]
    
    return nothing
end

@with_kw struct SchwarzschildChristoffelCache <: AbstractChristoffelCache
    l::Vector{Float64} = zeros(4)
    dH::Vector{Float64} = zeros(4)
    dl::Array{Float64, 2} = zeros(4,4)
    D::Array{Float64, 3} = zeros(4,4,4)
end

allocate_christoffel_cache(::SchwarzschildSpacetimeKerrSchildCoordinates) = SchwarzschildChristoffelCache()

function set_christoffel!(Γ, position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates, cache::SchwarzschildChristoffelCache)
    
    t, x, y, z = position
    M = spacetime.M

    r2 = x^2 + y^2 + z^2
    r  = sqrt(r2)
    r3 = r2*r

    #The scalar function and the null vector of the metric

    H = M/r

    l = cache.l
    l[1] = 1.
    l[2] = x/r
    l[3] = y/r 
    l[4] = z/r

    #Derivatives of the null vectors ( dl[a,b]=d_a l[b]). Indexes are down.

    dl = cache.dl
    fill!(dl, 0)

    dl[2,2] =1/r-x^2/r3
    dl[2,3] =-x*y/r3
    dl[2,4] =-x*z/r3

    dl[3,2] = dl[2,3]
    dl[3,3] = 1/r-y^2/r3
    dl[3,4] =-y*z/r3

    dl[4,2] = dl[2,4]
    dl[4,3] = dl[3,4]
    dl[4,4] = 1/r-z^2/r3

    #Derivatives of the scalar function H (dH[a]=d_a H). Index is down.

    dH = cache.dH
    fill!(dH, 0)

    dH[2] = -x*M/r3
    dH[3] = -y*M/r3
    dH[4] = -z*M/r3

    # Directional derivative of H in the direction of the null vector l  (l^a d_a H)
    l_dH = -M/r2

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


#Spherical coordinates

@with_kw struct SchwarzschildSpacetimeSphericalCoordinates <: AbstractSchwarzschildSpacetime 
    M::Float64
    @assert M >= 0.0
end

coordinates_topology(::SchwarzschildSpacetimeSphericalCoordinates) = SphericalTopology()
radius(position,::SchwarzschildSpacetimeSphericalCoordinates) = position[2]

function set_metric!(g, q, spacetime::SchwarzschildSpacetimeSphericalCoordinates)

    t, r, θ, φ = q

    M = spacetime.M

    fill!(g,0.0)
    g[1,1] = -(1-2M/r)
    g[2,2] =  1/(1-2M/r)
    g[3,3] =  r^2
    g[4,4] =  r^2*sin(θ)^2
    
    return nothing
end

function set_metric_inverse!(g, q, spacetime::SchwarzschildSpacetimeSphericalCoordinates)
    t, r, θ, φ = q
    M = spacetime.M
    fill!(g,0.0)
    g[1,1] = -1/(1-2M/r)
    g[2,2] =  1-2M/r
    g[3,3] =  1/r^2
    g[4,4] =  1/(r^2*sin(θ)^2)
    return nothing
end

allocate_christoffel_cache(::SchwarzschildSpacetimeSphericalCoordinates) = nothing

function set_christoffel!(Γ, point, spacetime::SchwarzschildSpacetimeSphericalCoordinates)
    t, r, θ, φ = point
    rs = 2*spacetime.M
    
    Γ[1,1,2] = rs/(2*r*(r-rs))          
    Γ[1,2,1] = Γ[1,1,2]    
    
    Γ[2,1,1] =  rs*(r-rs)/(2*r^3)
    Γ[2,2,2] = -rs/(2*r*(r-rs))
    Γ[2,3,3] = -(r-rs)
    Γ[2,4,4] = -(r-rs)*sin(θ)^2
    
    Γ[3,2,3] = 1.0/r   
    Γ[3,4,4] = -sin(θ)*cos(θ)
    Γ[3,3,2] = Γ[3,2,3]
    
    Γ[4,2,4] = 1.0/r
    Γ[4,3,4] = cot(θ)
    Γ[4,4,2] = Γ[4,2,4]
    Γ[4,4,3] = Γ[4,3,4]

    return nothing
end

#Common definitions

event_horizon_radius(spacetime::AbstractSchwarzschildSpacetime) = 2*spacetime.M
isco_radius(spacetime::AbstractSchwarzschildSpacetime) = 6*spacetime.M

function circular_geodesic_angular_speed(position, spacetime::AbstractSchwarzschildSpacetime, rotation_sense::AbstractRotationSense)
    M = spacetime.M
    r = radius(position, spacetime)
    Ω = sqrt(M)/r^1.5
    s = sign(rotation_sense)
    return s*Ω
end
