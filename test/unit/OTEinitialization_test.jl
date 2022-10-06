using Skylight, Test

@testset "initialization" begin
    
    @testset "configurations initialization" begin

        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)
    
        image_plane = ImagePlane(observer_distance = 1.0,
                                    observer_inclination_in_degrees = 137.0,
                                    horizontal_side_image_plane = 1.0,
                                    vertical_side_image_plane = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                   image_plane = image_plane,
                                                   initial_times = [0.0,1.0])
        
        rays = Skylight.my_zeros(configurations)
        @test sum(rays) == 0.0
        @test length(rays)/8 == 18
        @test Skylight.get_initial_times(configurations) == [0.0, 1.0]
    
    end

    @testset "four-momentum" begin

        @testset "set null ingoing past directed four-momentum" begin

            spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.99)

            position = [rand(),5.0,rand(),rand()]
            momentum = rand(4)
            momentum[1] = 0.0

            container = zeros(4,5)
            Skylight.dump_∂t_in!(container)
            Skylight.dump_metric_in!(container,position,spacetime)

            Skylight.set_null!(momentum,container)
            
            g = zeros(4,4)
            g = Skylight.set_metric!(g,position,spacetime)
            
            @test Skylight.norm_squared(momentum,g) ≈ 0.0 atol=1e-15
            
            @test momentum[1] == 1.0
            Skylight.set_ingoing_past_directed!(momentum)
            @test momentum[1] == -1.0
            
            momentum = rand(4)
            momentum[1] = 0.0

            Skylight.set_null_ingoing_past_directed!(momentum,container)
            
            @test Skylight.norm_squared(momentum,g) ≈ 0.0 atol=1e-15
            @test momentum[1] == -1.0
        end
    end

    @testset "single ray" begin
        
        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.0)

        image_plane = ImagePlane(observer_distance = sqrt(7.0),
                                    observer_inclination_in_degrees = 90.0,
                                    horizontal_side_image_plane = 1.0,
                                    vertical_side_image_plane = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                image_plane = image_plane,
                                                initial_times = [0.0,rand()])

        initial_time = configurations.initial_times[2]
        image_plane = configurations.image_plane
        coord_system = Skylight.coordinate_system_kind(configurations.spacetime)

        ray = zeros(8)
        pixel_coordinates = (1.0, 1.0)
        container = zeros(4,5)
        Skylight.dump_∂t_in!(container)

        Skylight.initialize_single!(ray,initial_time,pixel_coordinates,configurations,container)
        
        @test ray[1] ==  initial_time
        @test ray[2] ≈   sqrt(7.0)
        @test ray[3] ≈   1.0
        @test ray[4] ≈   1.0
        @test ray[5] == -1.0
        @test ray[6] ≈  -3*(-2*sqrt(7)+sqrt(69))/41
        @test ray[7] ≈   0.0
        @test ray[8] ≈   0.0  atol=1e-15
        
    end

    @testset "all rays" begin
        
        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.0)

        image_plane = ImagePlane(observer_distance = 3,
                                    observer_inclination_in_degrees = 90.0,
                                    horizontal_side_image_plane = 1.0,
                                    vertical_side_image_plane = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                image_plane = image_plane,
                                                initial_times = [0.1,1.5])

        rays = initialize(configurations)

        @views ray = rays[:,14]

        @test ray[1] ==  1.5
        @test ray[2] ≈   3.0
        @test ray[3] ≈   0.0
        @test ray[4] ≈   0.0 atol=1e-15
        @test ray[5] == -1.0
        @test ray[6] ≈  -0.2
        @test ray[7] ≈   0.0
        @test ray[8] ≈   0.0  atol=1e-15

        @views ray = rays[:,1]

        @test ray[1] ==  0.1
        @test ray[2] ≈   3.0
        @test ray[3] ≈  -0.5
        @test ray[4] ≈  -0.5 atol=1e-15
        @test ray[5] == -1.0
        @test ray[7] ≈   0.0
        @test ray[8] ≈   0.0  atol=1e-15
        
    end
end
