@testset "Cartesian coordinates" begin
    g = zeros(4, 4)
    spacetime = MinkowskiSpacetimeCartesianCoordinates()

    @test coordinates_topology(spacetime) == CartesianTopology()
    metric!(g, rand(4), spacetime)
    @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    cache = Skylight.ImagePlaneCache(spacetime)
    @views g1 = cache.metric
    metric!(g1, rand(4), spacetime)
    @test g1 == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    ginv = zeros(4, 4)
    metric_inverse!(ginv, rand(4), spacetime, g, nothing)
    @test g * ginv == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
end

@testset "Spherical coordinates" begin
    g = zeros(4, 4)
    spacetime = MinkowskiSpacetimeSphericalCoordinates()

    @test coordinates_topology(spacetime) == SphericalTopology()

    position = [rand(), 2.0, π / 2, rand()]
    metric!(g, position, spacetime)
    @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 4.0 0.0; 0.0 0.0 0.0 4.0]

    ginv = zeros(4, 4)
    metric_inverse!(ginv, position, spacetime, g, nothing)
    @test g * ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
end
