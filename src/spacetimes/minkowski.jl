"""
    MinkowskiSpacetimeCartesianCoordinates <: AbstractSpacetime

Minkowski spacetime using Cartesian coordinates. The metric takes the form:

See also: [Minkowski Spacetime](https://en.wikipedia.org/wiki/Minkowski_space). 
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

#Spherical coordinates
"""
    MinkowskiSpacetimeSphericalCoordinates <: AbstractSpacetime

Minkowski spacetime using spherical coordinates. The metric takes the form:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```
where `r` is the radial coordinate and `θ` (theta) is the polar angle.

See also: [Minkowski Spacetime](https://en.wikipedia.org/wiki/Minkowski_space). 
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
