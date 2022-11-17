using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

image_plane = ImagePlane(observer_distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side_image_plane = 20.0,
                         vertical_side_image_plane = 20.0,
                         horizontal_number_of_nodes = 300,
                         vertical_number_of_nodes = 300)

model = NovikovThorneDisk(inner_radius=6.0, outer_radius=18.0, rbound=2e-1)
        
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                   image_plane = image_plane,
                                   observed_times = [0.0],
                                   radiative_model = model)

solver_options = SolverOptions(method=VCABM(), reltol=1e-13, abstol=1e-21, output_type = SaveEndstate())

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations) #... or, define your own cb and cb_params

output_data = integrate_geodesics(initial_data, configurations, cb, cb_params, solver_options)