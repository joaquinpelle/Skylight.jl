@doc raw"""
    MinkowskiSpacetimeCartesianCoordinates <: AbstractSpacetime

[Minkowski Spacetime](https://en.wikipedia.org/wiki/Minkowski_space) in spherical coordinates. 

``ds^2 = -dt^2 + dx^2 + dy^2 + dz^2``

# Example
```
MinkowskiSpacetimeCartesianCoordinates()
```
"""
struct MinkowskiSpacetimeCartesianCoordinates <: AbstractSpacetime end

coordinates_topology(::MinkowskiSpacetimeCartesianCoordinates) = CartesianTopology()

function metric!(g, position, ::MinkowskiSpacetimeCartesianCoordinates)
    fill!(g, 0.0)
    g[1, 1] = -1.0
    g[2, 2] = 1.0
    g[3, 3] = 1.0
    g[4, 4] = 1.0
    return nothing
end

function metric_inverse!(g,
    q,
    spacetime::MinkowskiSpacetimeCartesianCoordinates,
    gaux,
    cache)
    metric!(g, q, spacetime)
end

function allocate_christoffel_cache(::MinkowskiSpacetimeCartesianCoordinates)
    return nothing
end

function christoffel!(Γ,
    position,
    ::MinkowskiSpacetimeCartesianCoordinates)
    fill!(Γ, 0.0)
    return nothing
end

#Spherical coordinates
@doc raw"""
    MinkowskiSpacetimeSphericalCoordinates <: AbstractSpacetime

[Minkowski Spacetime](https://en.wikipedia.org/wiki/Minkowski_space) in spherical coordinates. 

``ds^2 = -dt^2 + dr^2 + r^2 d\theta^2 + r^2 \sin^2 \theta d\phi^2``

where $r$ is the radial coordinate and $θ$ is the polar angle.

# Example
```
MinkowskiSpacetimeSphericalCoordinates()
```
"""
struct MinkowskiSpacetimeSphericalCoordinates <: AbstractSpacetime end

coordinates_topology(::MinkowskiSpacetimeSphericalCoordinates) = SphericalTopology()

function metric!(g, position, ::MinkowskiSpacetimeSphericalCoordinates)
    t, r, θ, φ = position
    fill!(g, 0.0)
    g[1, 1] = -1.0
    g[2, 2] = 1.0
    g[3, 3] = r^2
    g[4, 4] = r^2 * sin(θ)^2
    return nothing
end

function metric_inverse!(g, position, ::MinkowskiSpacetimeSphericalCoordinates, gaux, cache)
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

function christoffel!(Γ,
    position,
    ::MinkowskiSpacetimeSphericalCoordinates)
    fill!(Γ, 0.0)
    return nothing
end