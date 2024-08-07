using Skylight
using CairoMakie

spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.9)

camera = ImagePlane(distance = 500.0,
    observer_inclination_in_degrees = 90,
    horizontal_side = 23.0,
    vertical_side = 23.0,
    horizontal_number_of_pixels = 100,
    vertical_number_of_pixels = 100)

model = NovikovThorneDisk(inner_radius = 6.0, outer_radius = 18.0, Mdot_to_MEdd=0.1, η=0.1, M1=1.0)

configurations = VacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0)

initial_data = initialize(configurations)

cb, cbp = callback_setup(configurations; rhorizon_bound = 2e-1) #... or, define your own cb and cbp

run = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-6,
    abstol = 1e-6)

output_data = run.output_data

Iobs = observed_bolometric_intensities(initial_data, output_data, configurations)

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
