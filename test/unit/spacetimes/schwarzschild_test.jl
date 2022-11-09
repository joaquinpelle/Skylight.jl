@testset "Parameters" begin
        
    @test_throws AssertionError SchwarzschildSpacetimeKerrSchildCoordinates(M=-1.0)
    @test_throws AssertionError SchwarzschildSpacetimeSphericalCoordinates(M=-1.0)

end

@testset "Kerr-Schild coordinates" begin
        
    spacetime = SchwarzschildSpacetimeKerrSchildCoordinates(M=1.0)
    @test Skylight.coordinate_system_class(spacetime) == Skylight.CartesianClass()


    point = [rand(),0.0,3.0,4.0]

    g = zeros(4,4)
    Skylight.set_metric!(g,point,spacetime)
    
    @test g ≈ [-1.0+2/5 0.0 6/25 8/25; 0.0 1.0 0.0 0.0; 6/25 0.0 1.0+18/125 24/125; 8/25 0.0 24/125 1.0+32/125]
    
    ginv = zeros(4,4)
    Skylight.set_metric_inverse!(ginv,point,spacetime)

    @test ginv ≈ [-1.0-2/5 0.0 6/25 8/25; 0.0 1.0 0.0 0.0; 6/25 0.0 1.0-18/125 -24/125; 8/25 0.0 -24/125 1.0-32/125]
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
end

@testset "Spherical coordinates" begin

    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
    @test Skylight.coordinate_system_class(spacetime) == Skylight.SphericalClass()


    point = [rand(),5.0,π/3,0.0]

    g = zeros(4,4)
    Skylight.set_metric!(g,point,spacetime)
    
    @test g ≈ [-(1.0-2/5) 0.0 0.0 0.0; 0.0 1.0/(1.0-2/5) 0.0 0.0; 0.0 0.0 25.0 0.0; 0.0 0.0 0.0 25sin(π/3)^2]
    
    ginv = zeros(4,4)
    Skylight.set_metric_inverse!(ginv,point,spacetime)
    
    @test ginv ≈ [-1.0/(1.0-2/5) 0.0 0.0 0.0; 0.0 1.0-2/5 0.0 0.0; 0.0 0.0 1/25 0.0; 0.0 0.0 0.0 1/(25sin(π/3)^2)]
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

end

