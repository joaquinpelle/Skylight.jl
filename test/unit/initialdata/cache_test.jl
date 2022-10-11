using Skylight, Test

@testset "Cache" begin
    
    cache = Skylight.OTEInitialDataCache()

    vector = rand(4)
    Skylight.dump_vector_in!(cache, vector)
    @test cache.vector == vector
    @test cache.metric == zeros(4,4)

    Skylight.dump_∂t_in!(cache)
    @test cache.vector == [1.0, 0.0, 0.0, 0.0]
    @test cache.metric == zeros(4,4)

    spacetime = MinkowskiSpacetimeCartesianCoordinates()

    image_plane = ImagePlane(observer_distance = 1.0,
                                observer_inclination_in_degrees = 137.0,
                                horizontal_side_image_plane = 1.0,
                                vertical_side_image_plane = 1.0,
                                horizontal_number_of_nodes = 3,
                                vertical_number_of_nodes = 3)
    
    configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                            image_plane = image_plane,
                                            initial_times = [0.0,1.0])

    position = rand(4)
    Skylight.dump_metric_in!(cache,position,spacetime)
    
    @test cache.vector == [1.0, 0.0, 0.0, 0.0]
    @test cache.metric == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    metric, vector = Skylight.unpack_views(cache)

    @test metric == cache.metric
    @test vector == cache.vector

end