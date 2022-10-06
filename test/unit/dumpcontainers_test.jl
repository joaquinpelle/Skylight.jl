using Skylight, Test

@testset "Container" begin
    
    container = zeros(4,5)

    vector = rand(4)
    Skylight.dump_vector_in!(container, vector)
    @test container == [0 0 0 0 vector[1]; 0 0 0 0 vector[2]; 0 0 0 0 vector[3]; 0 0 0 0 vector[4]]

    Skylight.dump_∂t_in!(container)
    @test container == [0 0 0 0 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]

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
    Skylight.dump_metric_in!(container,position,spacetime)
    
    @test container == [-1.0 0.0 0.0 0.0 1.0; 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0]

    gμν, tμ = Skylight.unpack_views(container)

    @test gμν == container[:,1:4]
    @test tμ == container[:,5]

end