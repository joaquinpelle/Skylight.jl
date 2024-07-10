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

## Radiative model

```@docs
    rest_frame_four_velocity!(vector::AbstractVector, position::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, spacetime_cache::AbstractSpacetimeCache, model_cache::AbstractModelCache)
    rest_frame_bolometric_intensity(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)
    rest_frame_specific_intensity(position::AbstractVector, momentum::AbstractVector, energy::Real, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)
    is_final_position_at_source(position::AbstractVector, sapcetime::AbstractSpacetime, model::AbstractRadiativeModel)
    lorentz_factors(positions, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)
    allocate_cache(model::AbstractRadiativeModel)
```

### Surface emission models

```@docs
surface_differential!(differential::AbstractVector, position::AbstractVector, model::AbstractSurfaceEmissionModel, coords_top::AbstractCoordinatesTopology)
```

### Thermal emission models

```@docs
temperature(position::AbstractVector, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)
```

### Line emission models

```@docs
line_emission_profile(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)
```