using Skylight, Test

@testset "Synthetic polar cap" begin
    model = SyntheticPolarCap(star_radius = 5.0,
        angular_speed = 0.05,
        misalignment_angle_in_degrees = 90,
        angular_radius_in_degrees = 60,
        temperature = rand())

    @test model.misalignment_angle_in_radians ≈ π / 2
    @test model.angular_radius_in_radians ≈ π / 3

    spacetime = MinkowskiSpacetimeCartesianCoordinates()
    coords_top = CartesianTopology()
    points = Skylight.space_positions(10, spacetime, model, coords_top, nothing)

    for i in 1:10
        point = points[:, i]

        @test sum(point .* point) ≈ 25
        @test point[1] >= 5 * cos(π / 3)
    end

    vector = zeros(4)
    position = [rand(), 3.0, 0.0, 4.0]
    gμν = [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    model_cache = allocate_cache(model)
    rest_frame_four_velocity!(vector,
        position,
        gμν,
        spacetime,
        model,
        coords_top,
        model_cache)

    @test vector ≈ [1.0 / sqrt(0.9775), 0.0, 0.15 / sqrt(0.9775), 0.0]

    df = zeros(4)
    surface_differential!(df, position, model, coords_top)
    @test df == [0.0, 2 * position[2], 2 * position[3], 2 * position[4]]

    spacetime = KerrSpacetimeKerrSchildCoordinates(M = 1.0, a = 0.5)

    metric = zeros(4, 4)
    metric_inverse = zeros(4, 4)

    metric!(metric, position, spacetime)
    metric_inverse!(metric_inverse, position, spacetime, metric, nothing)

    normal = zeros(4)
    Skylight.unit_surface_normal!(normal,
        position,
        metric,
        metric_inverse,
        model,
        coords_top)

    tangent_vector = [0.0, -position[3], position[2], 0.0]

    @test Skylight.norm_squared(normal, metric) ≈ 1.0
    @test Skylight.scalar_product(normal, tangent_vector, metric)≈0.0 atol=1e-16

    vector = Skylight.lower_index(normal, metric)
    @test vector[2] / df[2] ≈ vector[4] / df[4]
end
