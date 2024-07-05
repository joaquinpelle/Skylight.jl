# Auxiliary functions

### Characteristic radii and geodesic angular speed

```@docs
radius(position::AbstractVector, spacetime::AbstractSpacetime)
event_horizon_radius(spacetime::AbstractBlackHoleSpacetime)
isco_radius(spacetime::AbstractSpacetime, rotation_sense::AbstractRotationSense)
mbco_radius(spacetime::AbstractSpacetime, rotation_sense::AbstractRotationSense)
circular_geodesic_angular_speed(position::AbstractVector, spacetime::AbstractSpacetime, rotation_sense::AbstractRotationSense)
```

### Spacetime symmetries

```@docs
stationarity(spacetime::AbstractSpacetime)
spherical_symmetry(spacetime::AbstractSpacetime)
axisymmetry(spacetime::AbstractSpacetime)
```

### Others

```@docs
equatorial_ring_areas(spacetime::AbstractSpacetime)
```