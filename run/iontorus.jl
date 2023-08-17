using Pkg
using Skylight
using CairoMakie

distance = 500
N = 400
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=0.5)
camera = PinholeCamera(position = [0.0, distance, π*(1/2-1/18), 0.0],
                        horizontal_aperture_in_degrees = rad2deg(atan(20/distance)),
                        vertical_aperture_in_degrees = rad2deg(atan(20/distance)),
                        horizontal_number_of_pixels = N,
                        vertical_number_of_pixels = N)
model = IonTorus(spacetime)
configurations = NonVacuumOTEConfigurations(spacetime=spacetime,
                                   camera = camera,
                                   radiative_model = model,
                                   unit_mass_in_solar_masses=1.0,
                                   observation_energies = exp10.(range(-10, stop=-5.5, length=20)))
initial_data = initialize(configurations)
cb, cbp = callback_setup(configurations; rhorizon_bound=2e-1) #... or, define your own cb and cbp
run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-8, abstol=1e-8)
output_data = run.output_data

#Image
xs, ys = axes_ranges(camera) 
zs = grid_view(output_data, configurations; energy_index = 1)

zs[zs.<1e-40] .= 1e-40
zs = log10.(zs)
logzmin = minimum(zs[zs.>-40])
logzmax = maximum(zs)

fig = Figure(font = "CMU Serif") #resolution=(600,400)
ax = Axis(fig[1,1], xlabel=L"\alpha", ylabel=L"\beta", ylabelsize = 26, xlabelsize = 26) 
hmap = heatmap!(xs, ys, zs; colormap=:gist_heat, interpolate=true, colorrange = (logzmax-10, logzmax))
Colorbar(fig[:, end+1], hmap, label=L"I", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
colsize!(fig.layout, 1, Aspect(1, 1.0))
colgap!(fig.layout, 7)
CairoMakie.save("torus_image.png", fig)

#Spectrum
F = spectrum(initial_data, output_data, configurations)
ν = erg_to_Hz(configurations.observation_energies)

fig = Figure()
ax = Axis(fig[1,1], xlabel=L"\nu \, [\text{Hz}]", ylabel=L"\nu F_{\nu} \,[\text{erg} \,\text{s}^{-1}\,\text{Hz}^{-1}]", xscale=log10, yscale=log10, xlabelsize = 26, ylabelsize = 26)
lines!(ax, ν, ν.*Hz_to_erg(F); linewidth=2.0, color=:blue)
CairoMakie.save("torus_spectrum.png", fig)
display(fig)

#Potential
potential = (x, z) -> Skylight.torus_normalized_potential([0.0,sqrt(x^2+z^2),acos(z/sqrt(x^2+z^2)),0.0], spacetime, model)

size = 10
N = 200
x_vals = LinRange(0.0, size, N)
z_vals = LinRange(0.0, size, N)
w_vals = [potential(x, z) for x in x_vals, z in z_vals]

w_vals[w_vals .< 0.0] .= 1e-40 
wmin = minimum(w_vals[w_vals .> 1e-40])
wmax = maximum(w_vals)

fig = Figure()
ax = Axis(fig[1, 1])
img = heatmap!(x_vals, z_vals, log10.(w_vals), colormap=cmap, interpolate=true, colorrange=(log10(wmin), log10(wmax)))

cbar = Colorbar(fig[1, 2], img, label="Normalized potential", ticklabelpad=10) # contour!(ax, z_vals, levels=10, linewidth=1) # Adding contour lines
CairoMakie.save("torus_potential.png", fig)
