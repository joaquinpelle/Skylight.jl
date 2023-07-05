using Skylight, Test

@testset "initialization" begin
    
    @testset "configurations initialization" begin

        spacetime = MinkowskiSpacetimeCartesianCoordinates()
    
        camera = ImagePlane(distance = 1.0,
                                    observer_inclination_in_degrees = 137.0,
                                    observation_times = [0.0,1.0],
                                    horizontal_side = 1.0,
                                    vertical_side = 1.0,
                                    horizontal_number_of_pixels = 3,
                                    vertical_number_of_pixels = 3)
        
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            radiative_model = DummyModel(),
                                                   camera = camera,
                                                   unit_mass_in_solar_masses=1.0)
        
        rays = Skylight.my_zeros(configurations)
        @test sum(rays) == 0.0
        @test length(rays)/8 == 18
    
    end

    @testset "four-momentum" begin

        @testset "set null ingoing past directed four-momentum" begin

            spacetime = MinkowskiSpacetimeCartesianCoordinates()
            position = [rand(),5.0,rand(),rand()]
            momentum = rand(4)
            momentum[1] = 0.0

            cache = Skylight.ImagePlaneCache(spacetime)
            Skylight.metric!(cache.metric,position,spacetime)
            Skylight.static_four_velocity!(cache)

            Skylight.set_null!(momentum,cache)
            
            g = zeros(4,4)
            metric!(g,position,spacetime)
            
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
        
        spacetime = MinkowskiSpacetimeCartesianCoordinates()

        camera = ImagePlane(distance = sqrt(7.0),
                                    observer_inclination_in_degrees = 90.0,
                                    observation_times = [0.0,rand()],
                                    horizontal_side = 3.0,
                                    vertical_side = 3.0,
                                    horizontal_number_of_pixels = 3,
                                    vertical_number_of_pixels = 3)
        
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            radiative_model = DummyModel(),
                                                camera = camera,
                                                unit_mass_in_solar_masses=1.0)

        initial_time = camera.observation_times[2]
        camera = configurations.camera
        coords_top = coordinates_topology(spacetime)
        

        ray = zeros(8)
        pixel_coordinates =  (1.0, 1.0)
        cache = Skylight.ImagePlaneCache(spacetime)
        cache.vector = [1.0, 0.0, 0.0, 0.0]

        Skylight.initialize_single!(ray,initial_time,pixel_coordinates,configurations,cache)
        
        @test ray[1] ==  initial_time
        @test ray[2] ≈   sqrt(7.0)
        @test ray[3] ≈   1.0
        @test ray[4] ≈   1.0
        @test ray[5] == -1.0
        @test ray[6] ≈  -1.0
        @test ray[7] ≈   0.0
        @test ray[8] ≈   0.0  atol=1e-15
        
    end

    @testset "all rays" begin
        
        spacetime = MinkowskiSpacetimeCartesianCoordinates()

        camera = ImagePlane(distance = 3,
                                    observer_inclination_in_degrees = 90.0,
                                    observation_times = [0.1,1.5],
                                    horizontal_side = 3.0,
                                    vertical_side = 3.0,
                                    horizontal_number_of_pixels = 3,
                                    vertical_number_of_pixels = 3)
        
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                            radiative_model = DummyModel(),
                                                camera = camera,
                                                unit_mass_in_solar_masses=1.0)

        rays = initialize(configurations)

        @views ray = rays[:,14]

        @test ray[1] ==  1.5
        @test ray[2] ≈   3.0
        @test ray[3] ≈   0.0
        @test ray[4] ≈   0.0 atol=1e-15
        @test ray[5] == -1.0
        @test ray[6] ≈  -1.0
        @test ray[7] ≈   0.0
        @test ray[8] ≈   0.0  atol=1e-15

        @views ray = rays[:,1]

        @test ray[1] ==  0.1
        @test ray[2] ≈   3.0
        @test ray[3] ≈  -1.0
        @test ray[4] ≈  -1.0 atol=1e-15
        @test ray[5] == -1.0
        @test ray[7] ≈   0.0
        @test ray[8] ≈   0.0  atol=1e-15
        
    end
end
