# Radiative models

Radiative models are represented by types that contain all the information necessary to calculate the radiation emitted by a source. The key functions defining these models include [`rest_frame_four_velocity!`](@ref), [`rest_frame_bolometric_intensity`](@ref), and [`rest_frame_specific_intensity`](@ref). These functions compute the four-velocity of the local rest frame, as well as the bolometric and specific intensity in that frame within the emission region.

## Examples

The following example demonstrates how to compute these quantities for a Novikov-Thorne accretion disk in Kerr spacetime at a given position. First, create a Kerr spacetime with mass `M = 1.0` and spin `a = 0.5`, and obtain its coordinates topology

```julia
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.5)
coords_top = coordinates_topology(spacetime)
```

Set the position at `t = 0.0`, `r = 3.0`, `θ = π/2` and `φ = 0.0`, and calculate the metric 

```julia
position = [0.0, 5.0, π/2, 0.0]  #t=0.0, r = 3.0, θ = π/2, φ = 0.0
g = zeros(4,4)
metric!(g, position, spacetime)
```

Create a Novikov-Thorne disk with the inner radius at the ISCO of the spacetime, and an outer radius of $1000$. Assume the unit mass is $10^{7}$ solar masses, and the accretion rate is $10\%$ of the Eddington accretion rate with a radiative efficiency of $10\%$:

```julia
disk = NovikovThorneDisk(inner_radius=isco_radius(spacetime, ProgradeRotation()), 
    outer_radius = 1000.0, 
    M1 = 1e7, 
    Mdot_to_MEdd = 0.1, 
    η = 0.1)
```

Compute the disk four-velocity at the position:

```julia
u = zeros(4)
rest_frame_four_velocity!(u, position, spacetime, disk, coords_top)
```

Since the emission in the local rest frame is isotropic, use a random photon momentum

```julia
momentum = rand(4) 
```

Compute the bolometric intensity and the specific intensity at an energy of $10^{-4}$ erg in the local rest frame:

```julia
Ibol = rest_frame_bolometric_intensity(position, momentum, u, g, spacetime, disk, coords_top)
energy = 1e-4
Ispec = rest_frame_specific_intensity(position, momentum, energy, u, g, spacetime, disk, coords_top)
```

For certain models, more specialized functions are available. For models with thermal emission, obtain the temperature at a given position with the [`temperature`](@ref) function:

```julia
T = temperature(position, spacetime, disk)
```

For line emission models, calculate the emissivity profile with the [`line_emission_profile`](@ref) function

```julia
disk = AccretionDiskWithFlatLamppostProfile(inner_radius=isco_radius(spacetime, ProgradeRotation()), 
    outer_radius = 1000.0,
    corona_height = 10.0)
ϵ = line_emission_profile(position, momentum, u, g, spacetime, disk, coords_top)
```

For models involving more complex calculations, cache objects are required to store intermediate results. Construct these caches with the  [`allocate_cache`](@ref) function.

```julia
spacetime_cache = allocate_cache(spacetime)
model_cache = allocate_cache(disk)
rest_frame_four_velocity!(u, position, spacetime, disk, coords_top, spacetime_cache, model_cache)
Ibol = rest_frame_bolometric_intensity(position, momentum, u, g, spacetime, disk, coords_top, model_cache)
Ispec = rest_frame_specific_intensity(position, momentum, energy, u, g, spacetime, disk, coords_top, model_cache)
```

## Catalogue

The currently available radiative models are as follows:

### Geometrically thin, optically thick accretion disks

These accretion disk models are hydrostationary, geometrically thin, and optically thick, assuming the disk is in the equatorial plane of the spacetime and occupying a finite radial range. The disk particles follow circular geodesics, requiring [`circular_geodesic_angular_speed`](@ref) to be implemented for the chosen spacetime.  

#### Thermal radiation

These models assume the disk emits as a blackbody with a position-dependent temperature.

```@docs
Skylight.ShakuraSunyaevDisk
```

```@docs
Skylight.NovikovThorneDisk
```

```@docs
Skylight.AccretionDiskWithTabulatedTemperature
```

```@docs
Skylight.RARDisk
```

#### Line emission

In these models, the radiation corresponds to line emission from the accretion disk.

```@docs
Skylight.AccretionDiskWithFlatLamppostProfile
```

```@docs
Skylight.AccretionDiskWithTabulatedProfile
```

### Geometrically thick, optically thin accretion disks

```@docs
Skylight.IonTorus
```

### Others

```@docs
Skylight.CircularHotSpot
```

```@docs
Skylight.LamppostCorona
```

```@docs
Skylight.StarAcrossWormhole
```