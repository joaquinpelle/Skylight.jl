using Skylight, Test

@testset "initialization" begin
    
    @testset "configurations initialization" begin

        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)
    
        image_plane = ImagePlane(distance = 1.0,
                                    observer_inclination_in_degrees = 137.0,
                                    horizontal_side = 1.0,
                                    vertical_side = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            radiative_model = DummyModel(),
                                                   image_plane = image_plane,
                                                   observed_times = [0.0,1.0],
                                                   unit_mass_in_solar_masses=1.0)
        
        rays = Skylight.my_zeros(configurations)
        @test sum(rays) == 0.0
        @test length(rays)/8 == 18
        @test Skylight.observed_times(configurations) == [0.0, 1.0]
    
    end

    @testset "four-momentum" begin

        @testset "set null ingoing past directed four-momentum" begin

            spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.99)

            position = [rand(),5.0,rand(),rand()]
            momentum = rand(4)
            momentum[1] = 0.0

            cache = Skylight.OTEInitialDataCache()
            Skylight.dump_∂t_in!(cache)
            Skylight.dump_metric_in!(cache,position,spacetime)

            Skylight.set_null!(momentum,cache)
            
            g = zeros(4,4)
            set_metric!(g,position,spacetime)
            
            @test Skylight.norm_squared(momentum,g) ≈ 0.0 atol=1e-15
            
            @test momentum[1] == 1.0
            Skylight.set_ingoing_past_directed!(momentum)
            @test momentum[1] == -1.0
            
            momentum = rand(4)
            momentum[1] = 0.0

            Skylight.set_null_ingoing_past_directed!(momentum,cache)
            
            @test Skylight.norm_squared(momentum,g) ≈ 0.0 atol=1e-15
            @test momentum[1] == -1.0
        end
    end

    @testset "single ray" begin
        
        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.0)

        image_plane = ImagePlane(distance = sqrt(7.0),
                                    observer_inclination_in_degrees = 90.0,
                                    horizontal_side = 1.0,
                                    vertical_side = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            radiative_model = DummyModel(),
                                                image_plane = image_plane,
                                                observed_times = [0.0,rand()],
                                                unit_mass_in_solar_masses=1.0)

        initial_time = configurations.observed_times[2]
        image_plane = configurations.image_plane
        coords_top = coordinates_topology(spacetime)
        

        ray = zeros(8)
        pixel_coordinates = (1.0, 1.0)
        cache = Skylight.OTEInitialDataCache()
        Skylight.dump_∂t_in!(cache)

        Skylight.initialize_single!(ray,initial_time,pixel_coordinates,configurations,cache)
        
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

        image_plane = ImagePlane(distance = 3,
                                    observer_inclination_in_degrees = 90.0,
                                    horizontal_side = 1.0,
                                    vertical_side = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            radiative_model = DummyModel(),
                                                image_plane = image_plane,
                                                observed_times = [0.1,1.5],
                                                unit_mass_in_solar_masses=1.0)

        rays = get_initial_data(configurations)

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
