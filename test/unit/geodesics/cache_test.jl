using Skylight, Test

@testset "Thread" begin
    
    spacetime = Skylight.KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)
    cache = Skylight.ThreadCache(christoffel_cache = Skylight.allocate_christoffel_cache(spacetime))

    @test cache.point == zeros(4)
    @test cache.velocity == zeros(4)
    @test cache.acceleration == zeros(4)
    @test cache.christoffel == zeros(4,4,4)
    @test cache.christoffel_cache.l == zeros(4)
    @test cache.christoffel_cache.dH == zeros(4)
    @test cache.christoffel_cache.dl == zeros(4,4)
    @test cache.christoffel_cache.D == zeros(4,4,4)


    cache = Skylight.allocate_multi_thread_cache(spacetime)
    @test length(cache) == Threads.nthreads()
    cache2 = cache[1]

    @test cache2.point == zeros(4)
    @test cache2.velocity == zeros(4)
    @test cache2.acceleration == zeros(4)
    @test cache2.christoffel == zeros(4,4,4)
    @test cache2.christoffel_cache.l == zeros(4)
    @test cache2.christoffel_cache.dH == zeros(4)
    @test cache2.christoffel_cache.dl == zeros(4,4)
    @test cache2.christoffel_cache.D == zeros(4,4,4)
    
end

@testset "Multi thread" begin
   
    spacetime = Skylight.KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)

    model = SyntheticPolarCap(number_of_points=10, 
                          star_radius=5.0,
                          angular_speed = 0.05, 
                          misalignment_angle_in_degrees=90,
                          angular_radius_in_degrees=60, 
                          temperature=1.0)

    image_plane = ImagePlane(observer_distance = 50.0,
                            observer_inclination_in_degrees = 45,
                            horizontal_side_image_plane = 30.0,
                            vertical_side_image_plane = 40.0,
                            horizontal_number_of_nodes = 50,
                            vertical_number_of_nodes = 50)   

    configurations = OTEConfigurations(spacetime=spacetime,
                                image_plane = image_plane,
                                initial_times = [0.0,1.0],
                                emission_model = model)
    rmax = 1.1*sqrt(6000)
    
    cb_params = Skylight.get_cb_params(model, configurations)
    
    geo_cache = Skylight.allocate_geodesics_cache(spacetime, cb_params)

    @test geo_cache.spacetime == spacetime
    @test geo_cache.cb_params == cb_params
    @test length(geo_cache.multi_thread) == Threads.nthreads()
    @test typeof(geo_cache.multi_thread[1]) == Skylight.ThreadCache{Skylight.KerrKSChristoffelCache}

end