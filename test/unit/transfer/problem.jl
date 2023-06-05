using Skylight, Test

@testset "Collect output" begin
    
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)

    camera = ImagePlane(distance = 500.0,
                            observer_inclination_in_degrees = 45,
                            horizontal_side = 10.0,
                            vertical_side = 10.0,
                            horizontal_number_of_pixels = 10,
                            vertical_number_of_pixels = 10)

    model = SyntheticPolarCap(
                            star_radius=5.0,
                            angular_speed = 0.05, 
                            misalignment_angle_in_degrees=90,
                            angular_radius_in_degrees=60, 
                            temperature=1.0)
            
    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                    camera = camera,
                                    observation_times = [0.0,1.0],
                                    radiative_model = model,
                                    unit_mass_in_solar_masses=1.0)

    initial_data = get_initial_data(configurations)

    cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp

    ensembleprob = Skylight.ensemble_problem(initial_data, configurations, cbp)
    
    #Also consider EnsembleSplitThreads() for multipixels and EnsembleGPUArray() for GPU

    stats = @timed sim = solve(ensembleprob, VCABM(), EnsembleThreads(); reltol=1e-14, abstol=1e-21, callback = cb, trajectories = 100)

    Skylight.collect_output(sim)
    
end
