@testset "Cartesian coordinates" begin
        
    g = zeros(4,4)
    spacetime = MinkowskiSpacetimeCartesianCoordinates()

    @test Skylight.coordinate_system_kind(spacetime) == Skylight.CartesianKind()
    Skylight.set_metric!(g,rand(4),spacetime)
    @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    cache = Skylight.OTEInitialDataCache()
    @views g1 = cache.metric
    Skylight.set_metric!(g1,rand(4),spacetime)
    @test g1 == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    ginv = zeros(4,4)
    Skylight.set_metric_inverse!(ginv,rand(4),spacetime)
    @test g*ginv == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    
end

@testset "Spherical coordinates" begin
    g = zeros(4,4)
    spacetime = MinkowskiSpacetimeSphericalCoordinates()

    @test Skylight.coordinate_system_kind(spacetime) == Skylight.SphericalKind()

    position = [rand(),2.0,π/2,rand()]
    Skylight.set_metric!(g,position,spacetime)
    @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 4.0 0.0; 0.0 0.0 0.0 4.0]

    ginv = zeros(4,4)
    Skylight.set_metric_inverse!(ginv,position,spacetime)
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

end