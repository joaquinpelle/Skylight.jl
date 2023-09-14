abstract type AbstractSchwarzschildSpacetime <: AbstractBlackHoleSpacetime end

@doc raw"""
    SchwarzschildSpacetimeKerrSchildCoordinates <: AbstractSchwarzschildSpacetime

[Schwarzschild spacetime](https://en.wikipedia.org/wiki/Schwarzschild_metric) in Kerr-Schild coordinates. The parameter $M$ is the mass. The metric is

``g_{\mu \nu} = \eta_{\mu \nu} + H l_{\mu} l_{\nu}``

where $\eta_{\mu \nu}$ is the flat metric, $H=2M/r$, and $l_{\mu}=(1,x,y,z)/r$.

# Constructor
```
SchwarzschildSpacetimeKerrSchildCoordinates(M=1.0)
```
"""
@with_kw struct SchwarzschildSpacetimeKerrSchildCoordinates <:
                AbstractSchwarzschildSpacetime
    M::Float64
    @assert M>=0.0 "M must be non-negative"
end

coordinates_topology(::SchwarzschildSpacetimeKerrSchildCoordinates) = CartesianTopology()

@doc raw"""
    radius(position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates)

    ``r = \sqrt{x^2 + y^2 + z^2}``
""" 
function radius(position, ::SchwarzschildSpacetimeKerrSchildCoordinates)
    sqrt(position[2]^2 + position[3]^2 + position[4]^2)
end

function metric!(g, position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates)
    M = spacetime.M
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    H = 2M / r
    l1 = 1.0
    l2 = x / r
    l3 = y / r
    l4 = z / r
    g[1, 1] = -1.0 + H * l1 * l1
    g[1, 2] = 0.0 + H * l1 * l2
    g[1, 3] = 0.0 + H * l1 * l3
    g[1, 4] = 0.0 + H * l1 * l4
    g[2, 1] = g[1, 2]
    g[2, 2] = 1.0 + H * l2 * l2
    g[2, 3] = 0.0 + H * l2 * l3
    g[2, 4] = 0.0 + H * l2 * l4
    g[3, 1] = g[1, 3]
    g[3, 2] = g[2, 3]
    g[3, 3] = 1.0 + H * l3 * l3
    g[3, 4] = 0.0 + H * l3 * l4
    g[4, 1] = g[1, 4]
    g[4, 2] = g[2, 4]
    g[4, 3] = g[3, 4]
    g[4, 4] = 1.0 + H * l4 * l4
    return nothing
end

function metric_inverse!(g::AbstractMatrix,
    position::AbstractVector,
    spacetime::SchwarzschildSpacetimeKerrSchildCoordinates,
    ::AbstractMatrix)
    M = spacetime.M
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    H = 2M / r

    l1 = -1.0
    l2 = x / r
    l3 = y / r
    l4 = z / r

    g[1, 1] = -1.0 - H
    g[1, 2] = 0.0 - H * l1 * l2
    g[1, 3] = 0.0 - H * l1 * l3
    g[1, 4] = 0.0 - H * l1 * l4
    g[2, 1] = g[1, 2]
    g[2, 2] = 1.0 - H * l2 * l2
    g[2, 3] = 0.0 - H * l2 * l3
    g[2, 4] = 0.0 - H * l2 * l4
    g[3, 1] = g[1, 3]
    g[3, 2] = g[2, 3]
    g[3, 3] = 1.0 - H * l3 * l3
    g[3, 4] = 0.0 - H * l3 * l4
    g[4, 1] = g[1, 4]
    g[4, 2] = g[2, 4]
    g[4, 3] = g[3, 4]
    g[4, 4] = 1.0 - H * l4 * l4
    return nothing
end

@with_kw struct SchwarzschildChristoffelCache <: AbstractChristoffelCache
    l::Vector{Float64} = zeros(4)
    dH::Vector{Float64} = zeros(4)
    dl::Array{Float64, 2} = zeros(4, 4)
    D::Array{Float64, 3} = zeros(4, 4, 4)
end

function allocate_christoffel_cache(::SchwarzschildSpacetimeKerrSchildCoordinates)
    SchwarzschildChristoffelCache()
end

function christoffel!(Γ,
    position,
    spacetime::SchwarzschildSpacetimeKerrSchildCoordinates,
    cache::SchwarzschildChristoffelCache)
    t, x, y, z = position
    M = spacetime.M
    r2 = x^2 + y^2 + z^2
    r = sqrt(r2)
    r3 = r2 * r

    #The scalar function and the null vector of the metric
    H = M / r
    l = cache.l
    l[1] = 1.0
    l[2] = x / r
    l[3] = y / r
    l[4] = z / r

    #Derivatives of the null vectors ( dl[a,b]=d_a l[b]). Indexes are down.
    dl = cache.dl
    fill!(dl, 0)
    dl[2, 2] = 1 / r - x^2 / r3
    dl[2, 3] = -x * y / r3
    dl[2, 4] = -x * z / r3
    dl[3, 2] = dl[2, 3]
    dl[3, 3] = 1 / r - y^2 / r3
    dl[3, 4] = -y * z / r3
    dl[4, 2] = dl[2, 4]
    dl[4, 3] = dl[3, 4]
    dl[4, 4] = 1 / r - z^2 / r3

    #Derivatives of the scalar function H (dH[a]=d_a H). Index is down.
    dH = cache.dH
    fill!(dH, 0)
    dH[2] = -x * M / r3
    dH[3] = -y * M / r3
    dH[4] = -z * M / r3

    # Directional derivative of H in the direction of the null vector l  (l^a d_a H)
    l_dH = -M / r2
    # Tensor product of the null vector derivatives with the null vector. dlablc[a,b,c]= dl[a,b]*l[c] 
    # Derivatives of the products H*la*lb:  D[a,b,c]= d_a (H*lb*lc) (The order of fors gives the order of indexes)
    # This computation is equivalent to D[a,b,c]=dH[a]*l[b]*l[c]+H*dl[a,b]*l[c]+H*dl[a,c]*l[b]
    D = cache.D
    for i in 1:4
        for j in 1:4
            for k in 1:4
                D[i, j, k] = dH[i] * l[j] * l[k] + H * dl[i, j] * l[k] + H * dl[i, k] * l[j]
            end
        end
    end
    #Christoffel symbols
    for i in 1:4
        sign = i == 1 ? -1 : 1
        for j in 1:4
            for k in 1:4
                Γ[i, j, k] = sign * (D[j, k, i] + D[k, j, i] - D[i, j, k] +
                              2 * H * l_dH * l[i] * l[j] * l[k])
            end
        end
    end
    return nothing
end

@doc raw"""
    SchwarzschildSpacetimeSphericalCoordinates <: AbstractSchwarzschildSpacetime

[Schwarzschild spacetime](https://en.wikipedia.org/wiki/Schwarzschild_metric) in spherical coordinates. The metric is

``ds^2 = -(1-2M/r) dt^2 + (1-2M/r)^{-1} dr^2 + r^2 d\theta^2 + r^2 \sin^2 \theta d\phi^2``

# Constructor
```
SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
```
"""
@with_kw struct SchwarzschildSpacetimeSphericalCoordinates <: AbstractSchwarzschildSpacetime
    M::Float64
    @assert M>=0.0 "M must be non-negative"
end

coordinates_topology(::SchwarzschildSpacetimeSphericalCoordinates) = SphericalTopology()
radius(position, ::SchwarzschildSpacetimeSphericalCoordinates) = position[2]

function metric!(g, q, spacetime::SchwarzschildSpacetimeSphericalCoordinates)
    t, r, θ, φ = q

    M = spacetime.M

    fill!(g, 0.0)
    g[1, 1] = -(1 - 2M / r)
    g[2, 2] = 1 / (1 - 2M / r)
    g[3, 3] = r^2
    g[4, 4] = r^2 * sin(θ)^2

    return nothing
end

function metric_inverse!(g::AbstractMatrix,
    position::AbstractVector,
    spacetime::SchwarzschildSpacetimeSphericalCoordinates,
    ::AbstractMatrix)
    t, r, θ, φ = position
    M = spacetime.M
    fill!(g, 0.0)
    g[1, 1] = -1 / (1 - 2M / r)
    g[2, 2] = 1 - 2M / r
    g[3, 3] = 1 / r^2
    g[4, 4] = 1 / (r^2 * sin(θ)^2)
    return nothing
end

allocate_christoffel_cache(::SchwarzschildSpacetimeSphericalCoordinates) = nothing

function christoffel!(Γ, position, spacetime::SchwarzschildSpacetimeSphericalCoordinates)
    t, r, θ, φ = position
    rs = 2 * spacetime.M

    Γ[1, 1, 2] = rs / (2 * r * (r - rs))
    Γ[1, 2, 1] = Γ[1, 1, 2]

    Γ[2, 1, 1] = rs * (r - rs) / (2 * r^3)
    Γ[2, 2, 2] = -rs / (2 * r * (r - rs))
    Γ[2, 3, 3] = -(r - rs)
    Γ[2, 4, 4] = -(r - rs) * sin(θ)^2

    Γ[3, 2, 3] = 1.0 / r
    Γ[3, 4, 4] = -sin(θ) * cos(θ)
    Γ[3, 3, 2] = Γ[3, 2, 3]

    Γ[4, 2, 4] = 1.0 / r
    Γ[4, 3, 4] = cot(θ)
    Γ[4, 4, 2] = Γ[4, 2, 4]
    Γ[4, 4, 3] = Γ[4, 3, 4]

    return nothing
end

#Common definitions

mass(spacetime::AbstractSchwarzschildSpacetime) = spacetime.M
event_horizon_radius(spacetime::AbstractSchwarzschildSpacetime) = 2 * spacetime.M
isco_radius(spacetime::AbstractSchwarzschildSpacetime, ::AbstractRotationSense) = 6 * spacetime.M

function circular_geodesic_angular_speed(position,
    spacetime::AbstractSchwarzschildSpacetime,
    rotation_sense::AbstractRotationSense)
    M = spacetime.M
    r = radius(position, spacetime)
    Ω = sqrt(M) / r^1.5
    s = sign(rotation_sense)
    return s * Ω
end