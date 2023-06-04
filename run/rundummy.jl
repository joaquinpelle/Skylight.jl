using Skylight
using CairoMakie

spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=0.9)

camera = ImagePlane(distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side = 20.0,
                         vertical_side = 20.0,
                         horizontal_number_of_nodes = 50,
                         vertical_number_of_nodes = 50)

model = Skylight.DummyExtendedRegion()

configurations = NonVacuumOTEConfigurations(spacetime=spacetime,
                                   camera = camera,
                                   radiative_model = model,
                                   observed_times = [0.0],
                                   observed_energies = [1.0],
                                   unit_mass_in_solar_masses=1.0)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations; rhorizon_bound = 0.3) #... or, define your own cb and cb_params

run = integrate(initial_data, configurations, cb, cb_params; τmax=2.0, method=VCABM(), reltol=1e-13, abstol=1e-21)
output_data = output_data(run)

xs, ys = get_pixel_coordinates_vectors(configurations) 
zs = view_as_grid(output_data, configurations; E_idx = 1)

fig = Figure(font = "CMU Serif") #resolution=(600,400)
ax = Axis(fig[1,1], xlabel=L"\alpha", ylabel=L"\beta", ylabelsize = 26, xlabelsize = 26) 
hmap = heatmap!(xs, ys, zs; colormap=:gist_heat, interpolate=true)
Colorbar(fig[:, end+1], hmap, label=L"I", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
CairoMakie.save("plot.png", fig)

