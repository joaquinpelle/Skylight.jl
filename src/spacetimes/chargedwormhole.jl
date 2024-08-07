"""
    AbstractChargedWormholeSpacetime <: AbstractSpacetime
"""
abstract type AbstractChargedWormholeSpacetime <: AbstractSpacetime end

stationarity(::AbstractChargedWormholeSpacetime) = IsStationary()
spherical_symmetry(::AbstractChargedWormholeSpacetime) = IsSphericallySymmetric()

allocate_christoffel_cache(::AbstractChargedWormholeSpacetime) = nothing

"""
    ChargedWormholeSpacetimeSphericalCoordinates <: AbstractChargedWormholeSpacetime

[Charged wormhole spacetime](https://arxiv.org/pdf/2210.14106.pdf) 
in usual spherical coordinates.

# Constructor
```
ChargedWormholeSpacetimeSphericalCoordinates(b0=1.0, Q=0.5)
```
"""
@with_kw struct ChargedWormholeSpacetimeSphericalCoordinates <:
                AbstractChargedWormholeSpacetime
    b0::Float64
    Q::Float64

    @assert b0 >= 0.0
    @assert abs(Q) < b0
end

coordinates_topology(::ChargedWormholeSpacetimeSphericalCoordinates) = SphericalTopology()

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::ChargedWormholeSpacetimeSphericalCoordinates)
    r = position[2]
    θ = position[3]

    b0 = spacetime.b0
    Q = spacetime.Q

    b = b0^2 / r

    gtt = -(1 + Q^2 / r^2)
    grr = 1 / (1 - b / r + Q^2 / r^2)

    g[1, 1] = gtt
    g[1, 2] = 0.0
    g[1, 3] = 0.0
    g[1, 4] = 0.0
    g[2, 1] = 0.0
    g[2, 2] = grr
    g[2, 3] = 0.0
    g[2, 4] = 0.0
    g[3, 1] = 0.0
    g[3, 2] = 0.0
    g[3, 3] = r^2
    g[3, 4] = 0.0
    g[4, 1] = 0.0
    g[4, 2] = 0.0
    g[4, 3] = 0.0
    g[4, 4] = r^2 * sin(θ)^2

    return nothing
end

function metric_inverse!(g::AbstractMatrix,
    position::AbstractVector,
    spacetime::ChargedWormholeSpacetimeSphericalCoordinates,
    ::AbstractMatrix)

    r = position[2]
    θ = position[3]

    b0 = spacetime.b0
    Q = spacetime.Q

    b = b0^2 / r

    gtt = -(1 + Q^2 / r^2)
    grr = 1 / (1 - b / r + Q^2 / r^2)

    g[1, 1] = 1.0 / gtt
    g[1, 2] = 0.0
    g[1, 3] = 0.0
    g[1, 4] = 0.0
    g[2, 1] = 0.0
    g[2, 2] = 1.0 / grr
    g[2, 3] = 0.0
    g[2, 4] = 0.0
    g[3, 1] = 0.0
    g[3, 2] = 0.0
    g[3, 3] = 1.0 / r^2
    g[3, 4] = 0.0
    g[4, 1] = 0.0
    g[4, 2] = 0.0
    g[4, 3] = 0.0
    g[4, 4] = 1.0 / (r^2 * sin(θ)^2)

    return nothing
end

function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::ChargedWormholeSpacetimeSphericalCoordinates)

    #Spacetime coordinates
    r = position[2]
    θ = position[3]

    b0 = spacetime.b0
    Q = spacetime.Q

    b = b0^2 / r
    db = -b0^2 / r^2

    gtt = -(1 + Q^2 / r^2)
    grr = 1 / (1 - b / r + Q^2 / r^2)

    dgtt = 2 * Q^2 / r^3
    dgrr = grr^2 * (db / r - b / r^2 + 2 * Q^2 / r^3)

    #Metric functions and derivatives
    dα = dgtt / (2 * gtt)
    dβ = dgrr / (2 * grr)

    #The Christoffel symbols
    Γ[1, 1, 2] = dα
    Γ[1, 2, 1] = Γ[1, 1, 2]

    Γ[2, 1, 1] = -gtt / grr * dα
    Γ[2, 2, 2] = dβ
    Γ[2, 3, 3] = -r / grr
    Γ[2, 4, 4] = -r / grr * sin(θ)^2

    Γ[3, 2, 3] = 1.0 / r
    Γ[3, 4, 4] = -sin(θ) * cos(θ)
    Γ[3, 3, 2] = Γ[3, 2, 3]

    Γ[4, 2, 4] = 1.0 / r
    Γ[4, 3, 4] = cos(θ) / sin(θ)
    Γ[4, 4, 2] = Γ[4, 2, 4]
    Γ[4, 4, 3] = Γ[4, 3, 4]

    return nothing
end

"""
    ChargedWormholeSpacetimeRegularCoordinates <: AbstractChargedWormholeSpacetime

[Charged wormhole spacetime](https://arxiv.org/pdf/2210.14106.pdf) 
in regular (throat crossing) coordinates.

# Constructor
```
ChargedWormholeSpacetimeRegularCoordinates(b0=1.0, Q=0.5)
```
"""
@with_kw struct ChargedWormholeSpacetimeRegularCoordinates <:
                AbstractChargedWormholeSpacetime
    b0::Float64
    Q::Float64

    @assert b0>=0.0 "b0 must be non-negative"
    @assert abs(Q)<b0 "|Q| must be smaller than b0"
end

function radius(position, spacetime::ChargedWormholeSpacetimeRegularCoordinates)
    t, l, θ, φ = position

    b0 = spacetime.b0
    Q = spacetime.Q
    return sqrt(l^2 + b0^2 - Q^2)
end

coordinates_topology(::ChargedWormholeSpacetimeRegularCoordinates) = SphericalTopology()

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::ChargedWormholeSpacetimeRegularCoordinates)
    t, l, θ, φ = position

    b0 = spacetime.b0
    Q = spacetime.Q

    r = sqrt(l^2 + b0^2 - Q^2)

    gtt = -(1 + Q^2 / r^2)
    gll = 1.0

    g[1, 1] = gtt
    g[1, 2] = 0.0
    g[1, 3] = 0.0
    g[1, 4] = 0.0
    g[2, 1] = 0.0
    g[2, 2] = gll
    g[2, 3] = 0.0
    g[2, 4] = 0.0
    g[3, 1] = 0.0
    g[3, 2] = 0.0
    g[3, 3] = r^2
    g[3, 4] = 0.0
    g[4, 1] = 0.0
    g[4, 2] = 0.0
    g[4, 3] = 0.0
    g[4, 4] = r^2 * sin(θ)^2

    return nothing
end

function metric_inverse!(g::AbstractMatrix,
    position::AbstractVector,
    spacetime::ChargedWormholeSpacetimeRegularCoordinates,
    ::AbstractMatrix)
    t, l, θ, φ = position

    b0 = spacetime.b0
    Q = spacetime.Q

    r = sqrt(l^2 + b0^2 - Q^2)

    gtt = -(1 + Q^2 / r^2)
    gll = 1.0

    g[1, 1] = 1.0 / gtt
    g[1, 2] = 0.0
    g[1, 3] = 0.0
    g[1, 4] = 0.0
    g[2, 1] = 0.0
    g[2, 2] = 1.0 / gll
    g[2, 3] = 0.0
    g[2, 4] = 0.0
    g[3, 1] = 0.0
    g[3, 2] = 0.0
    g[3, 3] = 1.0 / r^2
    g[3, 4] = 0.0
    g[4, 1] = 0.0
    g[4, 2] = 0.0
    g[4, 3] = 0.0
    g[4, 4] = 1.0 / (r^2 * sin(θ)^2)

    return nothing
end

function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::ChargedWormholeSpacetimeRegularCoordinates)

    #Spacetime coordinates
    t, l, θ, φ = position

    b0 = spacetime.b0
    Q = spacetime.Q

    r = sqrt(l^2 + b0^2 - Q^2)

    gtt = -(1 + Q^2 / r^2)
    gθθ = r^2
    gφφ = r^2 * sin(θ)^2

    dr_dl = l / r

    dlgtt = 2 * Q^2 / r^3 * dr_dl
    dlgθθ = 2 * l
    dlgφφ = 2 * l * sin(θ)^2
    dθgφφ = r^2 * 2 * sin(θ) * cos(θ)

    #The Christoffel symbols

    Γ[1, 1, 2] = 0.5 * dlgtt / gtt
    Γ[1, 2, 1] = Γ[1, 1, 2]

    Γ[2, 1, 1] = -dlgtt / 2
    Γ[2, 3, 3] = -dlgθθ / 2
    Γ[2, 4, 4] = -dlgφφ / 2

    Γ[3, 2, 3] = 0.5 * dlgθθ / gθθ
    Γ[3, 4, 4] = -sin(θ) * cos(θ)
    Γ[3, 3, 2] = Γ[3, 2, 3]

    Γ[4, 2, 4] = 0.5 * dlgφφ / gφφ
    Γ[4, 3, 4] = cos(θ) / sin(θ)
    Γ[4, 4, 2] = Γ[4, 2, 4]
    Γ[4, 4, 3] = Γ[4, 3, 4]

    return nothing
end
