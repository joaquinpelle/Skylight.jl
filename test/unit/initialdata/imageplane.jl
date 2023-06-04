using Skylight, Test

@testset "image plane" begin
    
    @testset "image plane area" begin

        camera = ImagePlane(distance = 1.0,
                                    observer_inclination_in_degrees = 137.0,
                                    horizontal_side = 1.0,
                                    vertical_side = 1.0,
                                    horizontal_number_of_pixels = 2,
                                    vertical_number_of_pixels = 2)

        @test Skylight.pixel_area(camera) == 0.25
        @test Skylight.area(camera) == 1.0

    end

    @testset "image plane" begin
        
        camera = ImagePlane(distance = 1.0,
                                    observer_inclination_in_degrees = 90.0,
                                    horizontal_side = 1.0,
                                    vertical_side = 1.0,
                                    horizontal_number_of_pixels = 2,
                                    vertical_number_of_pixels = 2)
        
        @test camera.observer_inclination_in_radians ≈ π/2
        @test Skylight.pixel_area(camera) == 0.25
        @test Skylight.area(camera) == 1.0

        iterator1 = Skylight.get_pixel_coordinates(camera)
        iterator2 = Iterators.product([-0.25, 0.25], [-0.25, 0.25])

        @test collect(iterator1) == collect(iterator2)

    end
end