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

    geo_cache = Skylight.allocate_geodesics_cache(spacetime)

    @test geo_cache.spacetime == spacetime
    @test length(geo_cache.multi_thread) == Threads.nthreads()
    @test typeof(geo_cache.multi_thread[1]) == Skylight.ThreadCache{Skylight.KerrKSChristoffelCache}

end