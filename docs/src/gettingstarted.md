# Getting started 

First, start a Julia REPL with `julia -t NT` where `NT` is the number of threads you want to use. Alternatively, you can copy the contents below into a script and run it using the same flag for multithreading.   

To get started, you will need a spacetime, a radiative model, and a camera. For example, to instantiate a Kerr spacetime in Kerr-Schild coordinates with mass $M=1$ and spin $a/M=0.5$, 

```julia
using Skylight

spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=0.5)
```

A catalogue of currently available spacetimes is at [Catalogue of spacetimes](@ref). Next, instantiate a radiative model as, e.g.

```julia
disk = NovikovThorneDisk(inner_radius=isco_radius(spacetime, ProgradeRotation()), outer_radius = 15.0)
```

See the currently available radiative models at [Catalogue of radiative models](@ref). Then, you can construct a camera as 

```julia
camera = PinholeCamera(position = [0.0, 500, π/2-π/20, 0.0],
                        horizontal_aperture_in_degrees = rad2deg(80/500),
                        vertical_aperture_in_degrees = rad2deg(80/500),
                        horizontal_number_of_pixels = 600,
                        vertical_number_of_pixels = 600)
```
See [Pinhole camera](@ref) for an explanation of this camera setup. Finally, gather these
objects into a configurations object. This is a vacuum transport problem, so use

```julia
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                        radiative_model=disk,
                                        camera = camera,
                                        unit_mass_in_solar_masses = 1e7)
```
where `unit_mass_in_solar_masses` is the unit mass in solar masses which determines fully the problem units together with $c=G=1$, and OTE stands for the observer-to-emitter scheme. This paticular configurations type will get the specialized methods for transport in vacuum. For more general non vacuum problems, use, e.g.

```julia
configurations = NonVacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1e7,
    observation_energies = exp10.(range(-10, stop = -5.5, length = 20)))
```
where `observation_energies` is a vector of the observation energies in CGS. 

Then, create the initial data as

```julia
initial_data = initialize(configurations)
```

The initial data is a matrix that has the initial conditions for each ray as columns, where the first four components are the spacetime coordinates, and the last four are the components of the initial four-momentum in the coordinate frame. 

Before running the ray-tracing, you need to specify a callback to be called at each step of the equations integration. For the default callback, use

```julia 
cb, cbp = callback_setup(configurations; rhorizon_bound=0.1)
```

Notice for this setup an extra parameter has to be specified, which determines the minimum radial distance to which rays can approach the event horizon before terminating the geodesic integration (this is because no future-directed rays can exit the event horizon, so, conversely, no past-directed rays can reach it). See the details for each callback at [Callbacks](@ref).

You can also define your own callbacks. For more details, see [Callbacks](@ref) and [Event Handling](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/). Finally, you can integrate the equations with

```julia
sim = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-8,
    abstol = 1e-8)
```
You can choose any of the available [solver methods](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) from [DifferentalEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/).

The output data can be obtained as

```julia
output_data = sim.output_data
```

This matrix contains the final coordinates and momenta of each ray, with the same structure as the initial data. Finally, you can compute, for instance, the observed bolometric intensity of the radiation field and produce an image as

```julia
using CairoMakie

Iobs = observed_bolometric_intensities(initial_data, output_data, configurations)

xs, ys = axes_ranges(camera)
zs = grid_view(Iobs, configurations)

fig = Figure(font = "CMU Serif")
ax = Axis(fig[1, 1],
    xlabel = L"\alpha \, [\mathrm{deg}]",
    ylabel = L"\beta \, [\mathrm{deg}]",
    ylabelsize = 26,
    xlabelsize = 26)
hmap = heatmap!(rad2deg.(xs), rad2deg.(ys), zs; colormap = :gist_heat, interpolate = true)
Colorbar(fig[:, end + 1],
    hmap,
    label = L"I \, [\mathrm{erg}/\mathrm{s}/\mathrm{cm^2}/\mathrm{sr}]",
    labelsize = 26,
    width    = 15,
    ticksize = 18,
    tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
display(fig)
```