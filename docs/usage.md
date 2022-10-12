Under construction...
#### Units

Skylight uses geometrized units $c = G = 1$.

#### Spacetimes

First, you have to choose a spacetime (and coordinate system). Currently, the available options are:

  * `MinkowskiSpacetimeCartesianCoordinates()`
  * `MinkowskiSpacetimeSphericalCoordinates()`
  * `SchwarzschildSpacetimeKerrSchildCoordinates()`
  * `SchwarzschildSpacetimeSphericalCoordinates()`
  * `KerrSpacetimeKerrSchildCoordinates()`
  * `KerrSpacetimeBoyerLindquistCoordinates()`
  * `JohannsenSpacetimeBoyerLindquistCoordinates()`

The following will be implemented soon:

  * `ChargedWormholeSphericalCoordinates()`
  * `ChargedWormholeRegularCoordinates()`

Some of the spacetimes require a set of parameters. For example, to construct a Kerr spacetime in Kerr-Schild coordinates with mass $M=1$ and spin $a/M=0.5$, run 

```
julia> spacetime = KerrSpacetimeKerrSchildCoordinate(M=1.0,a=0.5)
```

Spacetimes without parameters are constructed just like

```
julia> spacetime = MinkowskiSpacetimeCartesianCoordinates()
```

For more details about the spacetimes and the parameters they need, see the spacetimes documentation.

#### Emission models

The currently available emission model is:

  * `SyntheticPolarCap()`

The following will be implemented soon:

  * `BogdanovPolarCap()`
  * `OnionHotSpots()`
  * `ThinAccretionDisk()`
  * `BlackHoleCorona()`
  * `StarBehindWormhole()`

The emission models are determined by a set of parameters. To construct a synthetic polar cap, for example 

```
julia> model = Skylight.SyntheticPolarCap(number_of_points=10, 
                                          NS_radius=5.0,
                                          angular_speed = 0.05, 
                                          misalignment_angle_in_degrees=90,
                                          angular_radius_in_degrees=60, 
                                          temperature=1.0)
```

For the details, see the emission models documentation. 

#### Image plane

For the observet-to-emitter scheme, use the following to construct an image plane

```
julia> image_plane = ImagePlane(observer_distance = 500.0,
                                observer_inclination_in_degrees = 45,
                                horizontal_side_image_plane = 10.0,
                                vertical_side_image_plane = 10.0,
                                horizontal_number_of_nodes = 50,
                                vertical_number_of_nodes = 50)
```

#### Configurations

The initial data configurations in the observer-to-emitter scheme are constructed as follows

```
julia> configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                     image_plane = image_plane,
                                                     initial_times = [0.0,1.0])
```

In the emitter-to-observer scheme, use the following

```
julia> configurations = ETOInitialDataConfigurations(spacetime=spacetime,
                                                     emission_model=model,
                                                     number_of_packets_per_point = 100)
```

#### Initial data

Finally, for creating the initial data, use

```
julia> initial_data = initialize_data(configurations)
```