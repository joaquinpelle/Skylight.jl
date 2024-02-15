# Getting started 

### 1. Set up your environment


After you have installed the package, start a Julia REPL with multithreading optionally enabled as 

```bash
julia -t NT
``` 

where `NT` is the number of threads you want to use. Then, load the package as

```julia
using Skylight
```

### 2. Create a spacetime 

Here, we create a Kerr spacetime in Boyer-Lindquist coordinates: 

```julia
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=0.5)
```

Find other available spacetimes at [Catalogue of spacetimes](@ref). 

### 3. Create a radiative model

Here we create a Novikov-Thorne disk in prograde rotation around the black hole, specifying the inner and outer radii of the disk. The inner radius is set to the ISCO (Innermost Stable Circular Orbit) of the black hole, which depends on the spacetime parameters and rotation direction.

```julia
disk = NovikovThorneDisk(inner_radius=isco_radius(spacetime, ProgradeRotation()), outer_radius = 15.0)
```

Other available radiative models can be found at [Catalogue of radiative models](@ref). 

### 4. Set up a camera

Then, we create a pinhole camera, where the position is given by the spacetime coordinates of the observation point, the aperture is determined by the horizontal and vertical aperture in degrees, and the number of pixels in each direction are given. Here the apertures are set to capture a wide view of the accretion disk. For more details on the camera setup, see [Pinhole camera](@ref). 
```julia
camera = PinholeCamera(position = [0.0, 500, π/2-π/20, 0.0],
                        horizontal_aperture_in_degrees = rad2deg(80/500),
                        vertical_aperture_in_degrees = rad2deg(80/500),
                        horizontal_number_of_pixels = 600,
                        vertical_number_of_pixels = 600)
```

### 5. Set up the configurations object
We gather the spacetime, radiative model, and camera into a configurations object. This setup is for a vacuum transport problem (in the observer-to-emitter scheme).

```julia
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                        radiative_model=disk,
                                        camera = camera,
                                        unit_mass_in_solar_masses = 1e7)
```
where `unit_mass_in_solar_masses` is the unit mass in solar masses which, together with $c=G=1$, fully determines your unit system. This latter choice will affect the interpretation of the varius quantities set before as, e.g. the mass of the Kerr spacetime, which now is understood to correspond to $10^7$ solar masses.

For more general non-vacuum transfer problems, use 

```julia
configurations = NonVacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1e7,
    observation_energies = exp10.(range(-10, stop = -5.5, length = 20)))
```
where `observation_energies` has to be a vector of observation energies in CGS units. 

### 6. Generate the initial data
Create the initial data for the transport problem with

```julia
initial_data = initialize(configurations)
```

The initial data will be a matrix having the initial conditions for each observation direction as columns, with the first four components being the spacetime coordinates (the same as the camera position), and the last four components are the components of the initial momentum in the corresponding coordinate frame. 

### 7. Define callbacks
We define a callback to be called at each step of the equations integration. This is generally used to stop the integration under certain conditions as getting far away from the source, or intersecting the emitting surface. The default callbacks can be set up with

```julia 
cb, cbp = callback_setup(configurations; rhorizon_bound=0.1)
```

where `cb` is the Callback object, and `cbp` contains the callback parameters. In this particular setup, an extra parameter has can be passed as an argument (`rhorizon_bound`), which determines the minimum radial distance to which rays can approximate the event horizon before terminating the integration (this is because no future-directed rays can exit the event horizon, so, conversely, no past-directed rays can reach it). 

Additionally, custom callbacks can be defined. For more details, see [Callbacks](@ref) and [Event Handling](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/). 

### 8. Integrate the equations
Finally, we integrate the radiative transfer and geodesic equations, choosing a solver method and setting the relative and absolute tolerances. 

```julia
sim = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-8,
    abstol = 1e-8)
```

In certain setups, you may require lower errors for the integration to remain stable, especially when dealing with interpolated spacetimes. The `integrate` function is a wrapper for the DifferentialEquations.jl package's `solve` function. As such, you can pass any additional keyword arguments that are accepted by the `solve`. In particular, any of the available [solver methods](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) can be used. For more information, see the [DifferentialEquations.jl documentation](https://docs.sciml.ai/DiffEqDocs/stable/).

The output data can be extracted as

```julia
output_data = sim.output_data
```

This matrix contains the final coordinates and momenta of each ray, with the same structure as the initial data.

### 9. Visualize the results

Finally, we compute, for instance, the observed bolometric intensity of the radiation field and produce an image as

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