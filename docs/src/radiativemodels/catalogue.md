# Catalogue of radiative models

## Geometrically thin, optically thick accretion disks

These accretion disk models are hydrostationary, (infinitesimally) geometrically-thin and optically thick. The disk is assumed to be in the equatorial plane of the spacetime, occupying a certain radial range. The particles of the disk are assumed to follow circular geodesics (so the [`circular_geodesic_angular_speed`](@ref) function must be implemented for the chosen spacetime).  

### Thermal radiation

In these models, the disk is assumed to emit as a blackbody with a position-dependent effective temperature (given by the [`temperature`](@ref) function).

#### Shakura & Sunyaev accretion disk

```@docs
Skylight.ShakuraSunyaevDisk
```

#### Novikov & Thorne accretion disk

```@docs
Skylight.NovikovThorneDisk
```

#### Tabulated temperature 

```@docs
Skylight.AccretionDiskWithTabulatedTemperature
```

#### RAR extension of Shakura & Sunyaev accretion disk

```@docs
Skylight.RARDisk
```

### Line emission

In these models, the radiation corresponds to line emission from the disk. The emissivity profile is given by the [`line_emission_profile`](@ref) function.

### Flat spacetime lamppost corona line emission profile 

```@docs
Skylight.AccretionDiskWithFlatLamppostProfile
```

### Tabulated line emission profile 

```@docs
Skylight.AccretionDiskWithTabulatedProfile
```

## Geometrically thick, optically thin accretion disks

#### Ion torus 

```@docs
Skylight.IonTorus
```

## Others

### Circular hot spot on the surface of a neutron star

```@docs
Skylight.CircularHotSpot
```

### Lamppost corona

```@docs
Skylight.LamppostCorona
```

### Star across traversable wormhole

```@docs
Skylight.StarAcrossWormhole
```