# API

## Spacetime

### Geometric functions

#### Allocating methods
```@docs
metric(position::AbstractVector, spacetime::AbstractSpacetime)
metric_inverse(position::AbstractVector, spacetime::AbstractSpacetime)
volume_element(position::AbstractVector, spacetime::AbstractSpacetime)
christoffel(position::AbstractVector, spacetime::AbstractSpacetime)
```

#### Non-allocating methods

For computing these quantities within tight loops, we suggest to use the following non-allocating methods that write the output on preallocated arrays, and use caches to store intermediate results. 

```@docs
metric!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractSpacetimeCache})
metric_inverse!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
volume_element(position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
christoffel!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractChristoffelCache})
```

### Other functions

```@docs
coordinates_topology(spacetime::AbstractSpacetime)
radius(position::AbstractVector, spacetime::AbstractSpacetime)
event_horizon_radius(spacetime::AbstractBlackHoleSpacetime)
isco_radius(spacetime::AbstractSpacetime, rotation_sense::AbstractRotationSense)
mbco_radius(spacetime::AbstractSpacetime, rotation_sense::AbstractRotationSense)
circular_geodesic_angular_speed(position::AbstractVector, spacetime::AbstractSpacetime, rotation_sense::AbstractRotationSense)
```

### Cache allocation

```@docs
allocate_cache(spacetime::AbstractSpacetime)
allocate_christoffel_cache(spacetime::AbstractSpacetime)
```