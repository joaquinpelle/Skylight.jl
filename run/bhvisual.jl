using Skylight
using CairoMakie

spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.5)
model = NovikovThorneDisk(inner_radius = isco_radius(spacetime, ProgradeRotation()), outer_radius = 15.0)
distance = 200
cam = PinholeCamera(position = [0.0, distance, π / 2 - π / 20, 0.0],
    horizontal_aperture_in_degrees = 10,
    vertical_aperture_in_degrees = 10,
    horizontal_number_of_pixels = 1000,
    vertical_number_of_pixels = 1000)

configurations = VacuumOTEConfigurations(spacetime = spacetime,
    camera = cam,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0)

initial_data = initialize(configurations)

cb, cbp = callback_setup(configurations; rhorizon_bound = 2e-3) #... or, define your own cb and cbp

run = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-8,
    abstol = 1e-8)

output_data = run.output_data

Iobs = observed_bolometric_intensities(initial_data, output_data, configurations)

xs, ys = axes_ranges(cam)

zs = grid_view(Iobs, configurations)

fig = Figure(font = "CMU Serif")
ax = Axis(fig[1, 1])
hmap = heatmap!(xs, ys, zs / maximum(zs); colormap = :gist_heat, interpolate = true)
hidedecorations!(ax)
colsize!(fig.layout, 1, Aspect(1, 1.0))
display(fig)
CairoMakie.save("plot.png", fig)