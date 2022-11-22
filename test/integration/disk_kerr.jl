using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

image_plane = ImagePlane(observer_distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side_image_plane = 20.0,
                         vertical_side_image_plane = 20.0,
                         horizontal_number_of_nodes = 300,
                         vertical_number_of_nodes = 300)

model = NovikovThorneDisk(inner_radius=6.0, outer_radius=18.0)
        
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                   image_plane = image_plane,
                                   observed_times = [0.0],
                                   radiative_model = model)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations; rhorizon_bound=2e-1) #... or, define your own cb and cb_params

output_data = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)