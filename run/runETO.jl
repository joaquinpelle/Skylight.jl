using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

model = SyntheticPolarCap(star_radius=5.0,
                          angular_speed = 0.05, 
                          misalignment_angle_in_degrees=90,
                          angular_radius_in_degrees=60, 
                          temperature=1.0)

configurations = VacuumETOConfigurations(spacetime=spacetime,
                                   radiative_model = model,
                                   number_of_points=10,
                                   number_of_packets_per_point = 100, observer_distance = 500.0,
                                   unit_mass_in_solar_masses=1.0)

initial_data = initialize(configurations)

cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp

run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-13, abstol=1e-21)
output_data = output_data(run)
