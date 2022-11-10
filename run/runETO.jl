using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

model = SyntheticPolarCap(star_radius=5.0,
                          angular_speed = 0.05, 
                          misalignment_angle_in_degrees=90,
                          angular_radius_in_degrees=60, 
                          temperature=1.0)

configurations = ETOConfigurations(spacetime=spacetime,
                                   emission_model = model,
                                   number_of_points=10,
                                   number_of_packets_per_point = 100, observer_distance = 500.0)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

solver_options = SolverOptions(method=VCABM(), reltol=1e-13, abstol=1e-21, output_type = SaveEndpoint())

output_data = integrate_geodesics(initial_data, configurations, cb, cb_params, solver_options)