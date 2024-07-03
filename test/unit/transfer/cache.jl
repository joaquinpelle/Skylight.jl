using Skylight, Test

@testset "Thread" begin
    spacetime = KerrSpacetimeKerrSchildCoordinates(M = 1.0, a = 0.5)
    cache = Skylight.VacuumThreadCache(spacetime_cache = allocate_cache(spacetime),
        christoffel_cache = allocate_christoffel_cache(spacetime))

    @test cache.acceleration == zeros(4)
    @test cache.christoffel == zeros(4, 4, 4)
    @test cache.christoffel_cache.l == zeros(4)
    @test cache.christoffel_cache.dH == zeros(4)
    @test cache.christoffel_cache.dl == zeros(4, 4)
    @test cache.christoffel_cache.D == zeros(4, 4, 4)

    cache = Skylight.vacuum_multi_thread_cache(spacetime)
    @test length(cache) == Threads.nthreads()
    cache2 = cache[1]

    @test cache2.acceleration == zeros(4)
    @test cache2.christoffel == zeros(4, 4, 4)
    @test cache2.christoffel_cache.l == zeros(4)
    @test cache2.christoffel_cache.dH == zeros(4)
    @test cache2.christoffel_cache.dl == zeros(4, 4)
    @test cache2.christoffel_cache.D == zeros(4, 4, 4)
end

@testset "Multi thread" begin
    spacetime = KerrSpacetimeKerrSchildCoordinates(M = 1.0, a = 0.5)

    model = CircularHotSpot(
        star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
        spin_frequency_in_Hz = geometrized_to_CGS(0.05, Dimensions.frequency, M1 = 1.4),
        center_colatitude_in_degrees = 90.0,
        angular_radius_in_radians = deg2rad(60.0),
        M1 = 1.4,
        temperature_in_keV = 0.35)

    camera = ImagePlane(distance = 50.0,
        observer_inclination_in_degrees = 45,
        observation_times = [0.0, 1.0],
        horizontal_side = 30.0,
        vertical_side = 40.0,
        horizontal_number_of_pixels = 50,
        vertical_number_of_pixels = 50)

    configurations = VacuumOTEConfigurations(spacetime = spacetime,
        camera = camera,
        radiative_model = model,
        unit_mass_in_solar_masses = 1.0)
    rmax = 1.1 * sqrt(6000)

    cbp = Skylight.callback_parameters(spacetime, model, configurations)

    geo_cache = Skylight.transfer_cache(configurations, cbp)

    @test geo_cache.spacetime == spacetime
    @test geo_cache.cbp == cbp
    @test length(geo_cache.multi_thread) == Threads.nthreads()
    @test typeof(geo_cache.multi_thread[1]) ==
          Skylight.VacuumThreadCache{Nothing, Skylight.KerrChristoffelCache}
end
