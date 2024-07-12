# Novikov-Thorne disk around a Kerr black hole

In this tutorial, we will simulate the emission from a Novikov-Thorne disk around a Kerr black hole. Make sure you have installed Skylight.jl and CairoMakie.jl for visualizing the data.

#### 1. Import the necessary packages

```julia
using Skylight
using CairoMakie
```

#### 2. Define the spacetime

We define the spacetime as a Kerr black hole with mass `M = 1.0` and spin `a = 0.9` in Boyer-Lindquist coordinates.

```julia
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.9)
```

#### 3. Define the radiative model 

We create a Novikov-Thorne disk in prograde rotation around the black hole, specifying the inner and outer radii of the disk. The inner radius is set to the innermost stable circular orbit (ISCO) of the black hole, which depends on the spacetime parameters and rotation sense. The accretion rate is set to 10% of the Eddington accretion rate with a radiative efficiency of 10%. The unit mass is set to $10^7$ solar masses.

```julia
disk = NovikovThorneDisk(inner_radius=isco_radius(spacetime, ProgradeRotation()), 
    outer_radius = 1000.0, 
    M1 = 1e7, 
    Mdot_to_MEdd = 0.1, 
    η = 0.1,
    rotation_sense = ProgradeRotation())
```

#### 4. Set up the camera

```julia 
camera = PinholeCamera(position = [0.0, 500, π/2 - π/20, 0.0],
    horizontal_aperture_in_degrees = 4,
    vertical_aperture_in_degrees = 4,
    horizontal_number_of_pixels = 200,
    vertical_number_of_pixels = 200)
```

#### 5. Create the initial data

```julia
initial_data = initialize(configurations)
```

#### 6. Set up the callback

```julia
cb, cbp = callback_setup(configurations; rhorizon_bound = 2e-3) #... or, define your own cb and cbp
```

#### 7. Run the simulation

```julia
run = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-8,
    abstol = 1e-8)
```

#### 8. Extract the results

```julia
output_data = run.output_data
Iobs = observed_bolometric_intensities(initial_data, output_data, configurations)
```

#### 9. Visualize the results

```julia
xs, ys = axes_ranges(camera)
zs = grid_view(Iobs, configurations)

fig = Figure(font = "CMU Serif")
ax = Axis(fig[1, 1],
    xlabel = L"\alpha",
    ylabel = L"\beta",
    ylabelsize = 26,
    xlabelsize = 26)
hmap = heatmap!(xs, ys, zs / maximum(zs); colormap = :gist_heat, interpolate = true)
Colorbar(fig[:, end + 1],
    hmap,
    label = L"I",
    labelsize = 26,
    width = 15,
    ticksize = 18,
    tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
CairoMakie.save("plot.png", fig)
```
