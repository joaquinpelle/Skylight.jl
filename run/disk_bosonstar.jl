using Skylight
using CairoMakie

include("bosonstar_parameters.jl")

spacetime = BosonStarSpacetime(a=a_LBS1,b=b_LBS1)

camera = ImagePlane(distance = 500.0,
                         observer_inclination_in_degrees = 85,
                         horizontal_side = 90.0,
                         vertical_side = 90.0,
                         horizontal_number_of_pixels = 300,
                         vertical_number_of_pixels = 300)

model = AccretionDiskWithTabulatedTemperature(inner_radius=rin_LBS1, outer_radius=rout_LBS1, filename="TempLBS1.dat")
        
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                         camera = camera,
                                         observed_times = [0.0],
                                         radiative_model = model,
                                         unit_mass_in_solar_masses=1.0)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

run = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)

output_data = output_data(run)

bolometric_intensities = observed_bolometric_intensities(initial_data, output_data, configurations)

xs,ys = get_pixel_coordinates_vectors(configurations)

zs = view_as_grid(bolometric_intensities, configurations)

fig = Figure(font = "CMU Serif")
ax = Axis(fig[1,1], xlabel=L"\alpha", ylabel=L"\beta", ylabelsize = 26, xlabelsize = 26) 
hmap = heatmap!(xs, ys, zs; colormap=:gist_heat, interpolate=true)
Colorbar(fig[:, end+1], hmap, label=L"I", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
CairoMakie.save("plot.png", fig)