using Skylight
using CairoMakie

spacetime = SchwarzschildSpacetimeSphericalCoordinates(M = 1.0)

model = CircularHotSpot(
    star_radius_in_km = 12,
    spin_frequency_in_Hz = 200,
    center_colatitude_in_degrees = 90.0,
    angular_radius_in_radians = 0.01,
    M1 = 1.4,
    temperature_in_keV = 0.35)

configurations = VacuumETOConfigurations(spacetime = spacetime,
    radiative_model = model,
    number_of_points = 50,
    number_of_packets_per_point = 20000, 
    max_radius = 500.0,
    unit_mass_in_solar_masses = 1.4)

initial_data = initialize(configurations)

cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp

run = integrate(initial_data,
    configurations,
    cb,
    cbp;
    method = VCABM(),
    reltol = 1e-13,
    abstol = 1e-21)

output_data = run.output_data

observation_energies = [keV_to_erg(1.0)]
phase, skymap = specific_flux_skymap(initial_data, 
    output_data, 
    configurations, 
    observation_energies, 
    Nθ = 15, 
    Nϕ = 30)

skymap ./= maximum(skymap)
#TODO see extracting a desired angle
light_curve = skymap[1, 10, :]

light_curve = circshift(light_curve, 10)

fig = Figure(resolution = (800, 800))
ax = Axis(fig[1, 1])
lines!(ax, phase, light_curve)
display(fig)