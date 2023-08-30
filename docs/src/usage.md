# How to use

### Units

Skylight uses geometrized units $c = G = 1$.

### Spacetimes

First, you have to choose a spacetime (and coordinates topology). Currently, the available options are:

  * `MinkowskiSpacetimeCartesianCoordinates`
  * `MinkowskiSpacetimeSphericalCoordinates`
  * `SchwarzschildSpacetimeKerrSchildCoordinates`
  * `SchwarzschildSpacetimeSphericalCoordinates`
  * `KerrSpacetimeKerrSchildCoordinates`
  * `KerrSpacetimeBoyerLindquistCoordinates`
  * `FRKerrSpacetime`
  * `JohannsenSpacetime`
  * `RARSpacetime`
  * `BosonStarSpacetime`
  * `ChargedWormholeSphericalCoordinates`
  * `ChargedWormholeRegularCoordinates`

Some spacetimes depends on a set of parameters and others have no parameters. For example, to instantiate a Kerr spacetime in Kerr-Schild coordinates with mass $M=1$ and spin $a/M=0.5$, run 

```julia
spacetime = KerrSpacetimeKerrSchildCoordinate(M=1.0,a=0.5)
```

Spacetimes without parameters are instantiated just like

```julia
spacetime = MinkowskiSpacetimeCartesianCoordinates()
```

For more details about the spacetimes and the parameters they need, see the Spacetimes section of the documentation.

#### Radiative models

The currently available radiative models are:

  * `ShakuraSunyaevAccretionDisk` (provisory temperature profile)
  * `NovikovThorneAccretionDisk` (provisory temperature profile)
  * `RARAccretionDisk`
  * `AccretionDiskWithTabulatedTemperature`
  * `SyntheticPolarCap`

The following will be implemented soon:

  * `StarAcrossWormhole`
  * `BogdanovPolarCap`
  * `OnionHotSpots`
  * `BlackHoleCorona`

Radiative models also depend on a set of parameters. To construct a synthetic polar cap, for example 

```julia
model = Skylight.SyntheticPolarCap(star_radius=5.0,
                                          angular_speed = 0.05, 
                                          misalignment_angle_in_degrees=90,
                                          angular_radius_in_degrees=60, 
                                          temperature=1.0)
```

For the details, see the radiative models documentation. 

#### Image plane

For the observet-to-emitter scheme, use the following to construct an image plane

```julia
camera = ImagePlane(distance = 500.0,
                                observer_inclination_in_degrees = 45,
                                horizontal_side = 10.0,
                                vertical_side = 10.0,
                                horizontal_number_of_pixels = 600,
                                vertical_number_of_pixels = 600,
                                observation_times = [0.0,1.0])
```
#### Pinhole camera 

More generally, when asymptotic flatness isn't valid, or simply the observer is not located far away
enough, a pinhole camera can be used. It can be constructed, for instance, as
```julia
camera = PinholeCamera(position = [0.0, 500, π/2-π/20, 0.0],
                        horizontal_aperture_in_degrees = rad2deg(315/500),
                        vertical_aperture_in_degrees = rad2deg(315/500),
                        horizontal_number_of_pixels = 600,
                        vertical_number_of_pixels = 600)
```

#### Configurations

The initial data configurations in the observer-to-emitter scheme are constructed as follows

```julia
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                        radiative_model=model,
                                        camera = camera,
                                        unit_mass_in_solar_masses = 1.0)
```
where `unit_mass_in_solar_masses` is the unit mass in solar masses which determines fully the problem
units together with `c=G=1`.

In the emitter-to-observer scheme, use the following

```
julia> configurations = VacuumETOConfigurations(spacetime=spacetime,
                                                     radiative_model=model,
                                                     number_of_points = 100
                                                     number_of_packets_per_point = 100, max_radius = 500.0)
```

#### Initial data

Finally, for creating the initial data, use

```julia
initial_data = initialize(configurations)
```

#### Geodesics

Or you can specify another callback or method by...