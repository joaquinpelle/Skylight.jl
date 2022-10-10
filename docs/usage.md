Under construction...
#### Units

Skylight uses geometrized units $c = G = 1$.

#### The spacetime

First, you have to choose a spacetime (and coordinate system). Currently, the available options are:

  * `MinkowskiSpacetimeCartesianCoordinates()`
  * `MinkowskiSpacetimeSphericalCoordinates()`
  * `KerrSpacetimeKerrSchildCoordinates()`

The following are under testing:

  * `SchwarzschildSpacetimeKerrSchildCoordinates()`
  * `SchwarzschildSpacetimeSphericalCoordinates()`
  * `KerrSpacetimeBoyerLindquistCoordinates()`
  * `JohannsenSpacetimeBoyerLindquistCoordinates()`
  * `ChargedWormholeSphericalCoordinates()`
  * `ChargedWormholeRegularCoordinates()`

Some of the spacetimes require a set of parameters. For example, to create a Kerr spacetime in Kerr-Schild coordinates with $M=1$ and $a=0.5$, the call should be 

```
spacetime = KerrSpacetimeKerrSchildCoordinate(M=1.0,a=0.5)
```

This works analogously for the rest of the spacetimes. Spacetimes which do not have parameters can be created as

```
spacetime = MinkowskiSpacetimeCartesianCoordinates()
```

For more details about each spacetime and its parameters see the spacetimes documentation

#### Image plane

To build an image plane, use

```
image_plane = ImagePlane(observer_distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side_image_plane = 10.0,
                         vertical_side_image_plane = 10.0,
                         horizontal_number_of_nodes = 50,
                         vertical_number_of_nodes = 50)
```

#### Configurations

Before creating the initial data, you need to gather the information in a configurations data structure, which also
requires the specification of the initial times you desire to have initial data.

```
configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                              image_plane = image_plane,
                                              initial_times = [0.0,1.0])
```

#### Initial data

Finally, for creating the initial data, use

```
initial_data = initialize_data(configurations)
```