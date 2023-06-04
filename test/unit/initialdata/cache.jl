using Skylight, Test

@testset "Cache" begin
    
    @testset "Unpack views" begin
        
        cache = Skylight.ImagePlaneCache()
        spacetime = MinkowskiSpacetimeCartesianCoordinates()
        position = rand(4)

        Skylight.set_metric_and_four_velocity!(cache,position,spacetime)

        metric, vector = Skylight.unpack_views(cache)

        @test metric == cache.metric
        @test vector == cache.vector

        cache = Skylight.ETOInitialDataCache()
        spacetime = MinkowskiSpacetimeCartesianCoordinates()

        model = SyntheticPolarCap(
                                            star_radius=5.0,
                                            angular_speed = 0.05, 
                                            misalignment_angle_in_degrees=90,
                                            angular_radius_in_degrees=60, 
                                            temperature=rand())
        configurations = VacuumETOConfigurations(spacetime = spacetime, radiative_model = model, number_of_points=10, number_of_packets_per_point = 10, observer_distance = 500.0, unit_mass_in_solar_masses=1.0)

        Skylight.set_metric_and_tetrad!(cache, position, configurations)

        metric, metric_inverse, time_vector, triad = Skylight.unpack_views(cache)

        @test metric == cache.metric
        @test metric_inverse == cache.metric_inverse
        @test time_vector == cache.tetrad[:,1]
        @test triad == cache.tetrad[:,2:4]

    end

    @testset "Dumps" begin

        @testset "Observer to emitter" begin

            cache = Skylight.ImagePlaneCache()
            Skylight.dump_∂t_in!(cache)
            @test cache.vector == [1.0, 0.0, 0.0, 0.0]
            @test cache.metric == zeros(4,4)

            spacetime = MinkowskiSpacetimeCartesianCoordinates()

            Skylight.set_metric_and_four_velocity!(cache,position,spacetime)
            @test cache.metric == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
            
        end

        @testset "Emitter to observer" begin

            @testset "Generic model" begin

                cache = Skylight.ETOInitialDataCache()

                spacetime = MinkowskiSpacetimeCartesianCoordinates()
                model = DummyExtendedRegion()

                configurations = VacuumETOConfigurations(spacetime = spacetime, radiative_model = model, number_of_points=10, number_of_packets_per_point = 10, observer_distance = 500.0, unit_mass_in_solar_masses=1.0)

                packets = Skylight.my_zeros(configurations)
                cache = Skylight.ETOInitialDataCache()

                position = rand(4)

                Skylight.set_metric_and_tetrad!(cache, position, configurations)

                @views tetrad = cache.tetrad

                @test cache.metric == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
                @test tetrad[:,1] ≈ [1.0, 0.0, 0.0, 0.0]
                @test Skylight.scalar_product(tetrad[:,1], tetrad[:,2], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,1], tetrad[:,3], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,1], tetrad[:,4], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,2], tetrad[:,2], cache.metric) ≈ 1.0 atol = 1e-13 
                @test Skylight.scalar_product(tetrad[:,2], tetrad[:,3], cache.metric) ≈ 0.0 atol = 1e-13 
                @test Skylight.scalar_product(tetrad[:,2], tetrad[:,4], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,3], tetrad[:,3], cache.metric) ≈ 1.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,3], tetrad[:,4], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,4], tetrad[:,4], cache.metric) ≈ 1.0 atol = 1e-13
            
            end
            
            @testset "Surface model" begin
                
                cache = Skylight.ETOInitialDataCache()

                spacetime = MinkowskiSpacetimeCartesianCoordinates()
                model = SyntheticPolarCap( 
                                                    star_radius=5.0,
                                                    angular_speed = 0.05, 
                                                    misalignment_angle_in_degrees=90,
                                                    angular_radius_in_degrees=60, 
                                                    temperature=rand())
                configurations = VacuumETOConfigurations(spacetime = spacetime, radiative_model = model, number_of_points=10, number_of_packets_per_point = 10, observer_distance = 500.0, unit_mass_in_solar_masses=1.0)

                packets = Skylight.my_zeros(configurations)
                cache = Skylight.ETOInitialDataCache()

                position = [rand(), 3.0, 0.0, 4.0]

                Skylight.set_metric_and_tetrad!(cache, position, configurations)

                @views tetrad = cache.tetrad

                @test cache.metric == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
                @test cache.metric_inverse == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
                @test tetrad[:,1] ≈ [1.0/sqrt(0.9775), 0.0, 0.15/sqrt(0.9775), 0.0]
                @test tetrad[:,2] ≈ [0.0, 0.6, 0.0, 0.8]
                @test Skylight.scalar_product(tetrad[:,2], tetrad[:,3], cache.metric) ≈ 0.0 atol = 1e-13 
                @test Skylight.scalar_product(tetrad[:,2], tetrad[:,4], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,3], tetrad[:,3], cache.metric) ≈ 1.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,3], tetrad[:,4], cache.metric) ≈ 0.0 atol = 1e-13
                @test Skylight.scalar_product(tetrad[:,4], tetrad[:,4], cache.metric) ≈ 1.0 atol = 1e-13
            
            end
        end
    end
end