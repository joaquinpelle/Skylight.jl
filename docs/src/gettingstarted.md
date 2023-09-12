# Getting started 

To get started, you will need a spacetime, a radiative model, and a camera. For example, to instantiate a Kerr spacetime in Kerr-Schild coordinates with mass $M=1$ and spin $a/M=0.5$, 

```julia
using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinate(M=1.0,a=0.5)
```

A catalogue of currently available spacetimes is at [Catalogue of spacetimes](@ref). Next, instantiate a radiative model as, e.g.

```julia
disk = NovikovThorneAccretionDisk(inner_radius=isco_radius(spacetime),outer_radius=15.0)
```
See the currently available radiative models at [Catalogue of radiative models](@ref). Then, you can construct a camera as 

```julia
camera = PinholeCamera(position = [0.0, 500, π/2-π/20, 0.0],
                        horizontal_aperture_in_degrees = rad2deg(315/500),
                        vertical_aperture_in_degrees = rad2deg(315/500),
                        horizontal_number_of_pixels = 600,
                        vertical_number_of_pixels = 600)
```
See [Pinhole camera](@ref) for an explanation of this camera setup. Finally, gather these
objects into a configurations object as, e.g.

```julia
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                        radiative_model=model,
                                        camera = camera,
                                        unit_mass_in_solar_masses = 1.0)
```
where `unit_mass_in_solar_masses` is the unit mass in solar masses which determines fully the problem
units together with `c=G=1`. These are the configurations that will get the specialized routines for transport in vacuum. For more general non vacuum problems, use

```julia
configurations = NonVacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0,
    observation_energies = exp10.(range(-10, stop = -5.5, length = 20)))
```

From now on this quick start guide will focus on vacuum problems. For creating the initial data, use

```julia
initial_data = initialize(configurations)
```

Before running the ray-tracing, you need to specify a callback to be called at each step of the equations integration. For the default callback, use

```julia 
cb, cbp = callback_setup(configurations)
```

For more details on callbacks, see [Callbacks](@ref) and the [Event Handling](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) page of the DifferentialEquations.jl package documentaton. Finally, you can integrate the equations with

```julia
sim = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-8,
    abstol = 1e-8)
```
You can choose any of the available [solver methods](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) from DifferentalEquations.jl.

The output data can be obtained as

```julia
output_data = sim.output_data
```

Finally, you can compute, for instance, the observed bolometric intensity of the radiation field and produce an image as

```julia
using CairoMakie

Iobs = observed_bolometric_intensities(initial_data, output_data, configurations)

xs, ys = axes_ranges(camera)
zs = grid_view(Iobs, configurations)

fig = Figure(font = "CMU Serif")
ax = Axis(fig[1, 1],
    xlabel = L"\alpha",
    ylabel = L"\beta",
    ylabelsize = 26,
    xlabelsize = 26)
hmap = heatmap!(xs, ys, zs; colormap = :gist_heat, interpolate = true)
Colorbar(fig[:, end + 1],
    hmap,
    label = L"I",
    labelsize = 26,
    width = 15,
    ticksize = 18,
    tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
```