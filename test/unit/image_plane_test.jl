using Skylight, Test

@testset "image plane" begin
    
    @testset "image plane area" begin

        image_plane = ImagePlane(observer_distance = 1.0,
                                    observer_inclination_in_degrees = 137.0,
                                    horizontal_side_image_plane = 1.0,
                                    vertical_side_image_plane = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)

        @test Skylight.pixel_area(image_plane) == 0.25
        @test Skylight.area(image_plane) == 2.25

    end

    @testset "image plane" begin
        spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,a=0.5))
        
        image_plane = ImagePlane(observer_distance = 1.0,
                                    observer_inclination_in_degrees = 90.0,
                                    horizontal_side_image_plane = 1.0,
                                    vertical_side_image_plane = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                image_plane = image_plane,
                                                initial_times = [0.0])
        
        @test image_plane.observer_inclination_in_radians ≈ π/2
        @test Skylight.pixel_area(image_plane) == 0.25
        @test Skylight.area(image_plane) == 2.25

        iterator1 = Skylight.get_pixel_coordinates(configurations)
        iterator2 = Iterators.product([-0.5, 0.0, 0.5], [-0.5, 0.0, 0.5])

        @test collect(iterator1) == collect(iterator2)

    end
end