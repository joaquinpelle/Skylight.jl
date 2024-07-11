# Spacetimes

## Overview

The spacetimes are represented by types that contain the information needed to compute geodesics and other geometric properties. Most spacetimes have parameters, such as mass and spin, that are set at construction time. The core function defining the spacetimes is [`metric`](@ref), which returns the metric tensor at a given position in the spacetime. Additionally, [`metric_inverse`](@ref) returns the inverse metric tensor, [`volume_element`](@ref) returns the volume element, and [`christoffel`](@ref) returns the Christoffel symbols. Other functions are available to compute the radius of a position, the event horizon radius, the innermost stable circular orbit (ISCO) radius, the marginally bound circular orbit (MBCO) radius, and the circular geodesic angular speed, when applicable. See the [API](@ref) for more details.

## Example

The following code snippet demonstrates how to compute the metric tensor and the Christoffel symbols at a given position in the Schwarzschild spacetime in spherical coordinates:

```julia
using Skylight

# Define the Schwarzschild spacetime with mass M = 1
spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
position = [0.0, 4.0, π/2, 0.0]  #t=0.0, r = 4.0, θ = π/2, φ = 0.0
g = metric(position, spacetime)
Γ = christoffel(position, spacetime)
```

The example uses the [`metric`](@ref) and [`christoffel`](@ref) functions, which allocate the output arrays each time they are called. For performance-critical applications, especially within tight loops, use the non-allocating methods [`metric!`](@ref) and [`christoffel!`](@ref). These methods use preallocated arrays and auxiliary caches to store results and intermediate variables. The caches can be constructed with [`allocate_cache`](@ref) and [`allocate_christoffel_cache`](@ref), respectively.

```julia
g = zeros(4, 4)
cache = allocate_cache(spacetime)
metric!(g, position, spacetime, cache)
Γ = zeros(4, 4, 4)
christoffel_cache = allocate_christoffel_cache(spacetime)
christoffel!(Γ, position, spacetime, christoffel_cache)
```

For simpler spacetimes that do not require caches, the cache constructors return `nothing`, and the cache argument can be omitted in the non-allocating methods. These methods are more efficient than their allocating counterparts as they avoid cretaing an array for the output.

```julia
metric!(g, position, spacetime)
christoffel!(Γ, position, spacetime)
```

## Catalogue 

Below are the currently implemented spacetimes:

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

The spacetimes in `Skylight` are defined in various coordinate systems. Most coordinates used in practice have either Cartesian or spherical topology. For instance, Boyer-Lindquist coordinates for Kerr spacetime have spherical topology, whereas Kerr-Schild coordinates have Cartesian topology. The [`coordinates_topology`](@ref) function returns the topology of spacetime's coordinates. The following coordinate topologies are currently supported:

```@docs
Skylight.CartesianTopology
Skylight.SphericalTopology
```

## Abstract types

The spacetimes in `Skylight` are organized in a hierarchy of types leveraging Julia's type system and multiple dispatch to define common behaviors. All `Skylight` spacetimes are subtypes of the abstract type [`AbstractSpacetime`](@ref). For example, `AbstractBlackHoleSpacetime` is a subtype of `AbstractSpacetime` and is used to define spacetimes with an event horizon. This common feature among all black hole spacetimes is used, for example, to define a common callback for all of them which prevents geodesic integration getting too close to the event horizon, thus avoiding numerical instabilities. `AbstractKerrSpacetime` is a subtype of `AbstractBlackHoleSpacetime` and represents a Kerr spacetime, with various concrete Kerr spacetimes in different coordinate systems as subtypes. The following are the abstract spacetime types defined in `Skylight`: 

```@docs
Skylight.AbstractSpacetime
Skylight.AbstractBlackHoleSpacetime
Skylight.AbstractRegularCompactObjectSpacetime
Skylight.AbstractMinkowskiSpacetime
Skylight.AbstractSchwarzschildSpacetime
Skylight.AbstractKerrSpacetime
Skylight.AbstractChargedWormholeSpacetime
```
