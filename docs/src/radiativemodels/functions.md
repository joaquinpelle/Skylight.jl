# Functions

## Generic

### Rest-frame four velocities

```@docs
    rest_frame_four_velocity!(vector::AbstractVector, position::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, spacetime_cache::AbstractSpacetimeCache, model_cache::AbstractModelCache)
    lorentz_factors(positions::AbstractVector, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)
```

### Rest-frame intensities

```@docs
    rest_frame_four_velocity!(vector::AbstractVector, position::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, spacetime_cache::AbstractSpacetimeCache, model_cache::AbstractModelCache)
    rest_frame_bolometric_intensity(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)
    rest_frame_specific_intensity(position::AbstractVector, momentum::AbstractVector, energy::Real, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)
```

## For surface emission models

```@docs
surface_differential!(differential::AbstractVector, position::AbstractVector, model::AbstractSurfaceEmissionModel, coords_top::AbstractCoordinatesTopology)
```

## For thermal emission models

```@docs
temperature(position::AbstractVector, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)
```

## For line emission models

```@docs
line_emission_profile(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)
```

## Cache allocation

```@docs
allocate_cache(model::AbstractRadiativeModel)
```

## Miscellaneous

```@docs
is_final_position_at_source(position::AbstractVector, sapcetime::AbstractSpacetime, model::AbstractRadiativeModel)
```