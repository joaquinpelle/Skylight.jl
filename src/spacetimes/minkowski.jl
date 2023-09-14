abstract type AbstractMinkowskiSpacetime <: AbstractSpacetime end

@doc raw"""
    MinkowskiSpacetimeCartesianCoordinates <: AbstractMinkowskiSpacetime

[Minkowski Spacetime](https://en.wikipedia.org/wiki/Minkowski_space) in spherical coordinates. 

``ds^2 = -dt^2 + dx^2 + dy^2 + dz^2``

# Constructor
```
MinkowskiSpacetimeCartesianCoordinates()
```
"""
struct MinkowskiSpacetimeCartesianCoordinates <: AbstractMinkowskiSpacetime end

coordinates_topology(::MinkowskiSpacetimeCartesianCoordinates) = CartesianTopology()

function metric!(g::AbstractMatrix, position::AbstractVector, ::MinkowskiSpacetimeCartesianCoordinates)
    fill!(g, 0.0)
    g[1, 1] = -1.0
    g[2, 2] = 1.0
    g[3, 3] = 1.0
    g[4, 4] = 1.0
    return nothing
end

function metric_inverse!(g::AbstractMatrix,
    position::AbstractVector,
    spacetime::MinkowskiSpacetimeCartesianCoordinates,
    ::AbstractMatrix)
    metric!(g, position, spacetime)
end

function allocate_christoffel_cache(::MinkowskiSpacetimeCartesianCoordinates)
    return nothing
end

function christoffel!(Γ::AbstractArray,
    position::AbstractVector,
    ::MinkowskiSpacetimeCartesianCoordinates)
    fill!(Γ, 0.0)
    return nothing
end

#Spherical coordinates
@doc raw"""
    MinkowskiSpacetimeSphericalCoordinates <: AbstractMinkowskiSpacetime

[Minkowski Spacetime](https://en.wikipedia.org/wiki/Minkowski_space) in spherical coordinates. 

``ds^2 = -dt^2 + dr^2 + r^2 d\theta^2 + r^2 \sin^2 \theta d\phi^2``

where $r$ is the radial coordinate and $θ$ is the polar angle.

# Constructor
```
MinkowskiSpacetimeSphericalCoordinates()
```
"""
struct MinkowskiSpacetimeSphericalCoordinates <: AbstractMinkowskiSpacetime end

coordinates_topology(::MinkowskiSpacetimeSphericalCoordinates) = SphericalTopology()

function metric!(g::AbstractMatrix, position::AbstractVector, ::MinkowskiSpacetimeSphericalCoordinates)
    t, r, θ, φ = position
    fill!(g, 0.0)
    g[1, 1] = -1.0
    g[2, 2] = 1.0
    g[3, 3] = r^2
    g[4, 4] = r^2 * sin(θ)^2
    return nothing
end

function metric_inverse!(g::AbstractMatrix,
    position::AbstractVector,
    ::MinkowskiSpacetimeSphericalCoordinates,
    ::AbstractMatrix)
    t, r, θ, φ = position
    fill!(g, 0.0)
    g[1, 1] = -1.0
    g[2, 2] = 1.0
    g[3, 3] = 1 / r^2
    g[4, 4] = 1 / (r^2 * sin(θ)^2)
    return nothing
end

function allocate_christoffel_cache(::MinkowskiSpacetimeSphericalCoordinates)
    return nothing
end

function christoffel!(Γ::AbstractArray,
    position::AbstractVector,
    ::MinkowskiSpacetimeSphericalCoordinates)
    fill!(Γ, 0.0)
    return nothing
end