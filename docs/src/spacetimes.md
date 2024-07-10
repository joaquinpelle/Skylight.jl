# Spacetimes

## Overview

The spacetimes are represented by types that contain the information needed to compute geodesics and other properties. Most spacetimes have parameters (for example, their mass and spin) that are set at construction time. The most important function defining the spacetimes is [`metric!`](@ref), which returns the metric tensor at a given position in the spacetime. Additionally, [`metric_inverse!`](@ref) returns the inverse metric tensor, [`volume_element`](@ref) returns the volume element, and [`christoffel!`](@ref) returns the Christoffel symbols. Other functions are available to compute the radius of a position, the event horizon radius, the innermost stable circular orbit (ISCO) radius, the marginally bound circular orbit (MBCO) radius, and the circular geodesic angular speed, when applicable. See the [API](@ref) for more details.

As an example, the following code snippet shows how to compute the metric tensor and the Christoffel symbols at a given position in the Schwarzschild spacetime in spherical coordinates:

```julia
using Skylight

# Define the Schwarzschild spacetime with mass M = 1
spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
position = [0.0, 3.0, π/2, 0.0]  #t=0.0, r = 3.0, θ = π/2, φ = 0.0
g = metric(position, spacetime)
Γ = christoffel(position, spacetime)
```

The previous example uses the [`metric`](@ref) and [`christoffel`](@ref) functions, which allocate the output arrays each time they are called. However, when these functions need to be called within tight loops, we suggest using the non-allocating methods [`metric!`](@ref) and [`christoffel!`](@ref). These functions use preallocated arrays and auxiliary caches to store the results and intermediate variables. The caches can be constructed with the [`allocate_cache`](@ref) and [`allocate_christoffel_cache`](@ref) functions, respectively.

```julia
g = zeros(4, 4)
cache = allocate_cache(spacetime)
metric!(g, position, spacetime, cache)
Γ = zeros(4, 4, 4)
christoffel_cache = allocate_christoffel_cache(spacetime)
christoffel!(Γ, position, spacetime, christoffel_cache)
```

Some spacetimes are sufficiently simple that they do not require caches to store intermediate results, either in the calculation of the metric, the Christoffel symbols or both. In these cases, the cache constructors just return `nothing`, and the cache can be omitted as an argument in the non-allocating methods. These methods are still more efficient than the allocating ones, as they avoid allocating the output arrays.

```julia
metric!(g, position, spacetime)
christoffel!(Γ, position, spacetime)
```

## Catalogue 

The following are the currently implemented spacetimes:

### Minkowski spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractMinkowskiSpacetime && isconcretetype(t)
```

### Schwarzschild spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractSchwarzschildSpacetime && isconcretetype(t)
```

### Kerr spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractKerrSpacetime && isconcretetype(t)
```

### Johannsen spacetime

```@docs
Skylight.JohannsenSpacetime
```

### f(R)-Kerr spacetime

```@docs
Skylight.FRKerrSpacetime
```

### Charged wormhole spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractChargedWormholeSpacetime && isconcretetype(t)
```

### RAR spacetime

```@docs
Skylight.RARSpacetime
```

### Boson star spacetime

```@docs
Skylight.BosonStarSpacetime
```

## Coordinates topology

The spacetimes in `Skylight` are defined in different coordinate systems. The [`coordinates_topology`](@ref) function returns the topology of the coordinates of a given spacetime. The following are the currently implemented coordinate topologies:

```@docs
Skylight.CartesianTopology
Skylight.SphericalTopology
```

## Abstract types

The spacetimes in `Skylight` are organized in a hierarchy of types that serve to leverage Julia's type system and multiple dispatch to define common behaviors. All `Skylight` spacetimes are subtypes of the abstract type [`AbstractSpacetime`](@ref). For example, the abstract type `AbstractBlackHoleSpacetime` is a subtype of `AbstractSpacetime` and is used to define spacetimes that have an event horizon. This common feature among all black hole spacetimes is used, for example, to define a common callback for all of them which prevents geodesic integration too close to an event horizon, avoiding numerical instabilities. In turn, the abstract type `AbstractKerrSpacetime` represents a Kerr spacetime and has the various concrete Kerr spacetimes, in different coordinate systems, as subtypes. For reference, the following are the abstract spacetime types defined in `Skylight`: 

```@docs
Skylight.AbstractSpacetime
Skylight.AbstractBlackHoleSpacetime
Skylight.AbstractRegularCompactObjectSpacetime
Skylight.AbstractMinkowskiSpacetime
Skylight.AbstractSchwarzschildSpacetime
Skylight.AbstractKerrSpacetime
Skylight.AbstractChargedWormholeSpacetime
```
