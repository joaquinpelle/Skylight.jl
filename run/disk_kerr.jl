using Skylight

spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

camera = ImagePlane(distance = 500.0,
                         observer_inclination_in_degrees = 45,
                         horizontal_side = 20.0,
                         vertical_side = 20.0,
                         horizontal_number_of_nodes = 300,
                         vertical_number_of_nodes = 300)

model = NovikovThorneDisk(inner_radius=6.0, outer_radius=18.0)
        
configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                   camera = camera,
                                   observed_times = [0.0],
                                   radiative_model = model,
                                   unit_mass_in_solar_masses=1.0)

initial_data = get_initial_data(configurations)

cb, cb_params = get_callback_and_params(configurations; rhorizon_bound=2e-1) #... or, define your own cb and cb_params

run = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)

output_data = output_data(run)