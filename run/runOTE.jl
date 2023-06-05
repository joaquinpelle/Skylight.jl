using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

camera = ImagePlane(distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side = 10.0,
                         vertical_side = 10.0,
                         horizontal_number_of_pixels = 50,
                         vertical_number_of_pixels = 50)

model = SyntheticPolarCap(number_of_points=10, 
                          star_radius=5.0,
                          angular_speed = 0.05, 
                          misalignment_angle_in_degrees=90,
                          angular_radius_in_degrees=60, 
                          temperature=1.0)
        
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                   camera = camera,
                                   observed_times = [0.0,1.0],
                                   radiative_model = model,
                                   unit_mass_in_solar_masses=1.0)

initial_data = get_initial_data(configurations)

cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp

run = integrate(initial_data, configurations, cb, cbp; method=VCABM(), reltol=1e-13, abstol=1e-21)
output_data = output_data(run)
