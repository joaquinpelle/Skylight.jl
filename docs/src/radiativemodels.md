# Radiative models

## Overview

The radiative models are represented by types containing all the information necessary to calculate the radiation emitted by the source. The most important functions defining the radiative models are [`rest_frame_four_velocity!`](@ref), [`rest_frame_bolometric_intensity`](@ref), [`rest_frame_specific_intensity`](@ref). These functions are used to compute the four velocity of the local rest frame, and the bolometric and specific intensity in that frame, respectively.

The following example shows how to compute these quantities for a Novikov-Thorne accretion disk in the Kerr spacetime at a given position. First, create a Kerr spacetime with mass `M = 1.0` and spin `a = 0.5`, and obtain its coordinates topology

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

Create a Novikov-Thorne disk with inner radius at the ISCO of the spacetime, and outer radius at 1000.0, with the unit mass assumed to be $10^{7}$ solar masses, and the accretion rate set to 10% of the Eddington accretion rate with a radiative efficiency of 10%.

```julia
disk = NovikovThorneDisk(inner_radius=isco_radius(spacetime, ProgradeRotation()), 
    outer_radius = 1000.0, 
    M1 = 1e7, 
    Mdot_to_MEdd = 0.1, 
    η = 0.1)
```

Compute the disk four-velocity at the position

```julia
u = zeros(4)
rest_frame_four_velocity!(u, position, spacetime, disk, coords_top)
```

Since the emission in the local rest frame is isotropic, we can use a random photon momentum

```julia
momentum = rand(4) 
```

Compute the bolometric intensity and the specific intensity at an energy of $10^{-4}$ erg in the local rest frame

```julia
Ibol = rest_frame_bolometric_intensity(position, momentum, u, g, spacetime, disk, coords_top)
energy = 1e-4
Ispec = rest_frame_specific_intensity(position, momentum, energy, u, g, spacetime, disk, coords_top)
```

## Catalogue

### Geometrically thin, optically thick accretion disks

These accretion disk models are hydrostationary, (infinitesimally) geometrically-thin and optically thick. The disk is assumed to be in the equatorial plane of the spacetime, occupying a certain radial range. The particles of the disk are assumed to follow circular geodesics (so the [`circular_geodesic_angular_speed`](@ref) function must be implemented for the chosen spacetime).  

#### Thermal radiation

In these models, the disk is assumed to emit as a blackbody with a position-dependent temperature.

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