using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

model = Skylight.SyntheticPolarCap(number_of_points=10, 
                                   star_radius=5.0,
                                   angular_speed = 0.05, 
                                   misalignment_angle_in_degrees=90,
                                   angular_radius_in_degrees=60, 
                                   temperature=1.0)

configurations = ETOConfigurations(spacetime=spacetime,
                                              emission_model = model,
                                              number_of_packets_per_point = 100)

initial_data = initialize_data(configurations)