using Skylight, Test

@testset "Collect output" begin
    spacetime = KerrSpacetimeKerrSchildCoordinates(M = 1.0, a = 0.9)

    camera = ImagePlane(distance = 500.0,
        observer_inclination_in_degrees = 45,
        observation_times = [0.0, 1.0],
        horizontal_side = 10.0,
        vertical_side = 10.0,
        horizontal_number_of_pixels = 3,
        vertical_number_of_pixels = 3)

        model = CircularHotSpot(
            star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
            spin_frequency_in_Hz = geometrized_to_CGS(0.05, Dimensions.frequency, M1 = 1.4),
            center_colatitude_in_degrees = 90.0,
            angular_radius_in_radians = deg2rad(60.0),
            M1 = 1.4,
            temperature_in_keV = 0.35)

    configurations = VacuumOTEConfigurations(spacetime = spacetime,
        camera = camera,
        radiative_model = model,
        unit_mass_in_solar_masses = 1.0)

    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations) #... or, define your own cb and cbp
    ensembleprob = Skylight.ensemble_problem(initial_data, configurations, cbp)

    #Also consider EnsembleSplitThreads() for multinodes and EnsembleGPUArray() for GPU

    stats = @timed sim = solve(ensembleprob,
        VCABM(),
        EnsembleThreads();
        reltol = 1e-14,
        abstol = 1e-21,
        callback = cb,
        trajectories = 9)

    Skylight.collect_output(sim)
end
