using Skylight, Test

@testset "Thread" begin
    
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)
    cache = Skylight.VacuumThreadCache(spacetime_cache = allocate_cache(spacetime), christoffel_cache = allocate_christoffel_cache(spacetime))

    @test cache.acceleration == zeros(4)
    @test cache.christoffel == zeros(4,4,4)
    @test cache.christoffel_cache.l == zeros(4)
    @test cache.christoffel_cache.dH == zeros(4)
    @test cache.christoffel_cache.dl == zeros(4,4)
    @test cache.christoffel_cache.D == zeros(4,4,4)


    cache = Skylight.vacuum_multi_thread_cache(spacetime)
    @test length(cache) == Threads.nthreads()
    cache2 = cache[1]

    @test cache2.acceleration == zeros(4)
    @test cache2.christoffel == zeros(4,4,4)
    @test cache2.christoffel_cache.l == zeros(4)
    @test cache2.christoffel_cache.dH == zeros(4)
    @test cache2.christoffel_cache.dl == zeros(4,4)
    @test cache2.christoffel_cache.D == zeros(4,4,4)
    
end

@testset "Multi thread" begin
   
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)

    model = SyntheticPolarCap( 
                          star_radius=5.0,
                          angular_speed = 0.05, 
                          misalignment_angle_in_degrees=90,
                          angular_radius_in_degrees=60, 
                          temperature=1.0)

    camera = ImagePlane(distance = 50.0,
                            observer_inclination_in_degrees = 45,
                            observation_times = [0.0,1.0],
                            horizontal_side = 30.0,
                            vertical_side = 40.0,
                            horizontal_number_of_pixels = 50,
                            vertical_number_of_pixels = 50)   

    configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                camera = camera,
                                radiative_model = model,
                                unit_mass_in_solar_masses=1.0)
    rmax = 1.1*sqrt(6000)
    
    cbp = Skylight.callback_parameters(spacetime, model, configurations)
    
    geo_cache = Skylight.transfer_cache(configurations, cbp)

    @test geo_cache.spacetime == spacetime
    @test geo_cache.cbp == cbp
    @test length(geo_cache.multi_thread) == Threads.nthreads()
    @test typeof(geo_cache.multi_thread[1]) == Skylight.VacuumThreadCache{Nothing, Skylight.KerrChristoffelCache}

end