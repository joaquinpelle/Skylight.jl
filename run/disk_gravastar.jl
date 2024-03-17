using Skylight
using CairoMakie

distance = 750
Nres = 1000
apf = 50

spacetime = GravastarSpacetime(M = 1.0, R = 3.0, M_rho = 0.5)
camera = PinholeCamera(position = [0.0, distance, π / 2 - π / 20, 0.0],
    horizontal_aperture_in_degrees = rad2deg(apf / distance),
    vertical_aperture_in_degrees = rad2deg(apf / distance),
    horizontal_number_of_pixels = Nres,
    vertical_number_of_pixels = Nres)
model = NovikovThorneDisk(inner_radius = 6.0, outer_radius = 18.0)

configurations = VacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0)

initial_data = initialize(configurations)
cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp

run = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-8,
    abstol = 1e-8)

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
display(fig)
CairoMakie.save("gravastar.png", fig)
