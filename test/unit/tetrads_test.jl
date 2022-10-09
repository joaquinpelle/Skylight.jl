using Skylight, Test

@testset "Tetrads" begin
    
    @testset "Set random triad" begin
        
        triad = ones(4,3)
        Skylight.set_random_triad!(triad)
        @test sum(triad[1,:]) == 0.0

    end

    @testset "Orthonormalization" begin

        spacetime = Skylight.KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.9)
        triad = zeros(4,3)
        metric = zeros(4,4)
        position = [rand(), 3.0, 4.0, 1.0]
        Skylight.set_metric!(metric, position, spacetime)
        
        time_vector = [1.0, 0.0, 0.0, 0.0]
        Skylight.normalize_timelike!(time_vector, metric)
        
        Skylight.set_random_triad!(triad)

        Skylight.orthonormalize!(triad, time_vector, metric)

        scalar_product = Skylight.scalar_product

        @test scalar_product(time_vector, triad[:,1], metric) ≈ 0.0 atol=1e-13
        @test scalar_product(time_vector, triad[:,2], metric) ≈ 0.0 atol=1e-13
        @test scalar_product(time_vector, triad[:,3], metric) ≈ 0.0 atol=1e-13
        @test scalar_product(triad[:,1], triad[:,1], metric) ≈ 1.0 atol=1e-13
        @test scalar_product(triad[:,1], triad[:,2], metric) ≈ 0.0 atol=1e-13
        @test scalar_product(triad[:,1], triad[:,3], metric) ≈ 0.0 atol=1e-13
        @test scalar_product(triad[:,2], triad[:,2], metric) ≈ 1.0 atol=1e-13
        @test scalar_product(triad[:,2], triad[:,3], metric) ≈ 0.0 atol=1e-13
        @test scalar_product(triad[:,3], triad[:,3], metric) ≈ 1.0 atol=1e-13
    
    end

    @testset "Local triad" begin
        
        @testset "Random triad" begin

            triad = zeros(4,3)
            position = [rand(), 3.0, 4.0, 1.0]
            metric = zeros(4,4)
            metric_inverse = zeros(4,4)
            spacetime = Skylight.KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.9)

            Skylight.set_metric!(metric, position, spacetime)
            Skylight.set_metric_inverse!(metric_inverse, position, spacetime)

            time_vector = [1.0, 0.0, 0.0, 0.0]
            Skylight.normalize_timelike!(time_vector, metric)

            model = Skylight.DummyExtendedRegion()
            coord_system = Skylight.coordinate_system_kind(spacetime)
            

            Skylight.set_local_triad!(triad, position, time_vector, metric, metric_inverse, model, coord_system)

            scalar_product = Skylight.scalar_product

            @test scalar_product(time_vector, triad[:,1], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(time_vector, triad[:,2], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(time_vector, triad[:,3], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,1], triad[:,1], metric) ≈ 1.0 atol=1e-13
            @test scalar_product(triad[:,1], triad[:,2], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,1], triad[:,3], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,2], triad[:,2], metric) ≈ 1.0 atol=1e-13
            @test scalar_product(triad[:,2], triad[:,3], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,3], triad[:,3], metric) ≈ 1.0 atol=1e-13
        
        end
        @testset "Surface emission triad" begin

            triad = zeros(4,3)
            position = [rand(), 3.0, 4.0, 0.0]
            metric = zeros(4,4)
            metric_inverse = zeros(4,4)
            spacetime = Skylight.KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.9)

            Skylight.set_metric!(metric, position, spacetime)
            Skylight.set_metric_inverse!(metric_inverse, position, spacetime)

            time_vector = [1.0, 0.0, 0.0, 0.0]
            Skylight.normalize_timelike!(time_vector, metric)

            model = Skylight.SyntheticPolarCap(number_of_points=10, 
                                    NS_radius=5.0,
                                    angular_speed = 0.05, 
                                    misalignment_angle_in_degrees=90,
                                    angular_radius_in_degrees=60, 
                                    temperature=rand())

            coord_system = Skylight.coordinate_system_kind(spacetime)
            
            Skylight.set_local_triad!(triad, position, time_vector, metric, metric_inverse, model, coord_system)

            scalar_product = Skylight.scalar_product

            @test scalar_product(time_vector, triad[:,1], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(time_vector, triad[:,2], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(time_vector, triad[:,3], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,1], triad[:,1], metric) ≈ 1.0 atol=1e-13
            @test scalar_product(triad[:,1], triad[:,2], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,1], triad[:,3], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,2], triad[:,2], metric) ≈ 1.0 atol=1e-13
            @test scalar_product(triad[:,2], triad[:,3], metric) ≈ 0.0 atol=1e-13
            @test scalar_product(triad[:,3], triad[:,3], metric) ≈ 1.0 atol=1e-13
        
        end
    end

    @testset "Component conversion" begin

        spacetime = Skylight.MinkowskiSpacetimeCartesianCoordinates()
        metric = zeros(4,4)
        position = rand(4)
        Skylight.set_metric!(metric, position, spacetime)
    
        tetrad = zeros(4,4)
        tetrad[:,1] .= [sqrt(4/3), 1/sqrt(3), 0.0, 0.0]
        tetrad[:,2] .= [1/sqrt(3), sqrt(4/3), 0.0, 0.0]
        tetrad[:,3] .= [0.0, 0.0, 1.0, 1.0]
        tetrad[:,4] .= [0.0, 0.0, 1.0,-1.0]
    
        @views begin
            time_vector = tetrad[:,1] 
            triad = tetrad[:,2:4]
        end
    
        Skylight.orthonormalize!(triad,time_vector,metric)
    
        momenta = zeros(4,4)
    
        momenta[:,1] = [1.0, 1.0, 0.0, 0.0]
        momenta[:,2] = [1.0, 2.0, 0.0, 1.0]
        momenta[:,3] = [0.0, 1.0, 1.0, 1.0]
        momenta[:,4] = [0.0, 0.0,-1.0, 3.0]
    
        Skylight.set_coordinate_components_from_tetrad_components!(momenta, tetrad)
    
        @test momenta[:,1] ≈ [sqrt(4/3)+1/sqrt(3), sqrt(4/3)+1/sqrt(3), 0.0, 0.0]
        @test momenta[:,2] ≈ [sqrt(4/3)+2/sqrt(3), 2*sqrt(4/3)+1/sqrt(3), 1/sqrt(2), -1/sqrt(2)]
        @test momenta[:,3] ≈ [1/sqrt(3),sqrt(4/3), 2/sqrt(2), 0.0]
        @test momenta[:,4] ≈ [0.0, 0.0, 2/sqrt(2), -4/sqrt(2)]
    
    end
end
