using Test
using Skylight

@testset "Tetrad components" begin
    @testset "Set time and triad components" begin
        Nvectors = 8
        kμ = zeros(4, Nvectors)

        Skylight.unit_time_components!(kμ)

        for i in 1:Nvectors
            @test kμ[1, i] ≈ 1.0
            @test sum(kμ[2:4, i] .* kμ[2:4, i]) == 0.0
        end

        model = DummyModel()
        trait = opaque_interior_surface_trait(model)
        Skylight.packets_unit_random_triad_components!(kμ, trait)

        for i in 1:Nvectors
            @test kμ[1, i] ≈ 1.0
            @test sum(kμ[2:4, i] .* kμ[2:4, i]) ≈ 1.0
        end

        model = CircularHotSpot(
            star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
            spin_frequency_in_Hz = geometrized_to_CGS(0.05/(2π), Dimensions.frequency, M1 = 1.4),
            center_colatitude_in_degrees = 90.0,
            angular_radius_in_radians = deg2rad(60.0),
            M1 = 1.4,
            temperature_in_keV = 0.35)

        trait = opaque_interior_surface_trait(model)
        Skylight.packets_unit_random_triad_components!(kμ, trait)

        for i in 1:Nvectors
            @test kμ[1, i] ≈ 1.0
            @test sum(kμ[2:4, i] .* kμ[2:4, i]) ≈ 1.0
            @test kμ[2, i] >= 0.0
        end
    end
end

@testset "Set position and momenta" begin
    spacetime = MinkowskiSpacetimeCartesianCoordinates()
    model = CircularHotSpot(
        star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
        spin_frequency_in_Hz = geometrized_to_CGS(0.05/(2π), Dimensions.frequency, M1 = 1.4),
        center_colatitude_in_degrees = 90.0,
        angular_radius_in_radians = deg2rad(60.0),
        M1 = 1.4,
        temperature_in_keV = 0.35)
    configurations = VacuumETOConfigurations(spacetime = spacetime,
        radiative_model = model,
        number_of_points = 10,
        number_of_packets_per_point = 10,
        max_radius = 500.0,
        unit_mass_in_solar_masses = 1.0)

    packets = Skylight.my_zeros(configurations)
    cache = Skylight.ETOInitialDataCache(spacetime, model)

    position = [rand(), 3.0, 0.0, 4.0]

    Skylight.metric_and_tetrad!(cache, position, configurations)

    @views tetrad = cache.tetrad

    @test cache.metric ==
          [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    @test cache.metric_inverse ==
          [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    @test tetrad[:, 1] ≈ [1.0 / sqrt(0.9775), 0.0, 0.15 / sqrt(0.9775), 0.0]
    @test tetrad[:, 2] ≈ [0.0, 0.6, 0.0, 0.8]
    @test Skylight.scalar_product(tetrad[:, 2], tetrad[:, 3], cache.metric)≈0.0 atol=1e-12
    @test Skylight.scalar_product(tetrad[:, 2], tetrad[:, 4], cache.metric)≈0.0 atol=1e-12
    @test Skylight.scalar_product(tetrad[:, 3], tetrad[:, 3], cache.metric)≈1.0 atol=1e-12
    @test Skylight.scalar_product(tetrad[:, 3], tetrad[:, 4], cache.metric)≈0.0 atol=1e-12
    @test Skylight.scalar_product(tetrad[:, 4], tetrad[:, 4], cache.metric)≈1.0 atol=1e-12

    @views xμ = packets[1:4, 1:10]
    Skylight.packets_positions!(xμ, position)

    for i in 1:10
        @test xμ[:, i] == position
    end

    @views kμ = packets[1:4, 1:10]

    Skylight.packets_momenta!(kμ, tetrad, model)

    for i in 1:10
        @test Skylight.norm_squared(kμ[:, i], cache.metric)≈0.0 atol=1e-12
        @test Skylight.scalar_product(kμ[:, i], cache.tetrad[:, 2], cache.metric) >= 0.0
    end
end

@testset "Get initial positions" begin
    spacetime = MinkowskiSpacetimeCartesianCoordinates()
    model = CircularHotSpot(
        star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
        spin_frequency_in_Hz = geometrized_to_CGS(0.05/(2π), Dimensions.frequency, M1 = 1.4),
        center_colatitude_in_degrees = 90.0,
        angular_radius_in_radians = deg2rad(60.0),
        M1 = 1.4,
        temperature_in_keV = 0.35)
    configurations = VacuumETOConfigurations(spacetime = spacetime,
        radiative_model = model,
        number_of_points = 3,
        number_of_packets_per_point = 3,
        max_radius = 500.0,
        unit_mass_in_solar_masses = 1.0)

    metric = zeros(4, 4)
    metric!(metric, zeros(4), spacetime)

    times = Skylight.zero_times(configurations)

    @test times == [0.0, 0.0, 0.0]

    positions = zeros(4, 3)

    @views space_pos = positions[2:4, :]

    space_pos .= Skylight.space_positions(configurations)

    for i in 1:3
        position = positions[:, i]
        @test position[1] == 0.0
        @test position[2]^2 + position[3]^2 + position[4]^2 ≈ 25
        @test Skylight.scalar_product(position, [0.0, 1.0, 0.0, 0.0], metric) >=
              5 * cos(π / 3)
    end
end

@testset "Packets at position" begin
    spacetime = MinkowskiSpacetimeCartesianCoordinates()
    model = CircularHotSpot(
        star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
        spin_frequency_in_Hz = geometrized_to_CGS(0.05/(2π), Dimensions.frequency, M1 = 1.4),
        center_colatitude_in_degrees = 90.0,
        angular_radius_in_radians = deg2rad(60.0),
        M1 = 1.4,
        temperature_in_keV = 0.35)

    configurations = VacuumETOConfigurations(spacetime = spacetime,
        radiative_model = model,
        number_of_points = 10,
        number_of_packets_per_point = 10,
        max_radius = 500.0,
        unit_mass_in_solar_masses = 1.0)

    packets = Skylight.my_zeros(configurations)
    cache = Skylight.initial_data_cache(configurations)
    position = [rand(), 3.0, 0.0, 4.0]

    @views packets_at_position = packets[:, 1:10]

    Skylight.initialize_packets_at_position!(packets_at_position,
        position,
        cache,
        configurations)

    @views begin
        xμ = packets_at_position[1:4, :]
        kμ = packets_at_position[5:8, :]
    end

    for i in 1:10
        @test xμ[:, i] == position
    end

    for i in 1:10
        @test Skylight.norm_squared(kμ[:, i], cache.metric)≈0.0 atol=1e-12
        @test Skylight.scalar_product(kμ[:, i], cache.tetrad[:, 2], cache.metric) >= 0.0
    end
end

@testset "Initialization" begin
    spacetime = MinkowskiSpacetimeCartesianCoordinates()
    model = CircularHotSpot(
        star_radius_in_km = 1e-5*geometrized_to_CGS(5.0, Dimensions.length, M1 = 1.4),
        spin_frequency_in_Hz = geometrized_to_CGS(0.05/(2π), Dimensions.frequency, M1 = 1.4),
        center_colatitude_in_degrees = 90.0,
        angular_radius_in_radians = deg2rad(60.0),
        M1 = 1.4,
        temperature_in_keV = 0.35)
    configurations = VacuumETOConfigurations(spacetime = spacetime,
        radiative_model = model,
        number_of_points = 4,
        number_of_packets_per_point = 3,
        max_radius = 500.0,
        unit_mass_in_solar_masses = 1.0)

    packets = initialize(configurations; tasks_per_thread = 1)

    metric = zeros(4, 4)
    metric!(metric, zeros(4), spacetime)

    for i in 1:4
        @views begin
            packets_at_position = packets[:, i:(i + 2)]
            xμ = packets_at_position[1:4, :]
            kμ = packets_at_position[5:8, :]
        end

        normal = zeros(4)

        for j in 1:3
            normal = xμ[:, j]
            Skylight.normalize_spacelike!(normal, metric)

            @test xμ[1, j] == 0.0
            @test Skylight.scalar_product(xμ[:, j], [0.0, 1.0, 0.0, 0.0], metric) >=
                  5 * cos(π / 3)
            @test Skylight.norm_squared(kμ[:, j], metric)≈0.0 atol=1e-12
            @test Skylight.scalar_product(kμ[:, j], normal, metric) >= 0.0
        end
    end
end
