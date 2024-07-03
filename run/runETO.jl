using Skylight

spacetime = SchwarzschildSpacetimeSphericalCoordinates(M = 1.0)

model = CircularHotSpot(
    star_radius_in_km = 12,
    angular_speed_in_Hz = 200,
    center_colatitude_in_degrees = 90.0,
    angular_radius_in_radians = 1,
    M1 = 1.4,
    temperature_in_keV = 0.35)

configurations = VacuumETOConfigurations(spacetime = spacetime,
    radiative_model = model,
    number_of_points = 10,
    number_of_packets_per_point = 100, max_radius = 500.0,
    unit_mass_in_solar_masses = 1.0)

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