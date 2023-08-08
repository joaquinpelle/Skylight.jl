@testset "Parameters" begin
        
    @test_throws AssertionError SchwarzschildSpacetimeKerrSchildCoordinates(M=-1.0)
    @test_throws AssertionError SchwarzschildSpacetimeSphericalCoordinates(M=-1.0)

end

@testset "Kerr-Schild coordinates" begin
        
    spacetime = SchwarzschildSpacetimeKerrSchildCoordinates(M=1.0)
    @test coordinates_topology(spacetime) == CartesianTopology()


    point = [rand(),0.0,3.0,4.0]

    g = zeros(4,4)
    metric!(g,point,spacetime)
    
    @test g ≈ [-1.0+2/5 0.0 6/25 8/25; 0.0 1.0 0.0 0.0; 6/25 0.0 1.0+18/125 24/125; 8/25 0.0 24/125 1.0+32/125]
    
    ginv = zeros(4,4)
    metric_inverse!(ginv,point,spacetime, g, nothing)

    @test ginv ≈ [-1.0-2/5 0.0 6/25 8/25; 0.0 1.0 0.0 0.0; 6/25 0.0 1.0-18/125 -24/125; 8/25 0.0 -24/125 1.0-32/125]
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    cache = allocate_christoffel_cache(spacetime)
    
    kerr_spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.0)
    kerr_cache = allocate_christoffel_cache(kerr_spacetime)
    
    position = rand(4)
    Γ = zeros(4,4,4)
    Γ_kerr = zeros(4,4,4)

    christoffel!(Γ,position,spacetime,cache)
    christoffel!(Γ_kerr,position,kerr_spacetime,kerr_cache)
    
    @test Γ ≈ Γ_kerr
end

@testset "Spherical coordinates" begin

    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
    @test coordinates_topology(spacetime) == SphericalTopology()


    point = [rand(),5.0,π/3,0.0]

    g = zeros(4,4)
    metric!(g,point,spacetime)
    
    @test g ≈ [-(1.0-2/5) 0.0 0.0 0.0; 0.0 1.0/(1.0-2/5) 0.0 0.0; 0.0 0.0 25.0 0.0; 0.0 0.0 0.0 25sin(π/3)^2]
    
    ginv = zeros(4,4)
    metric_inverse!(ginv,point,spacetime, g, nothing)
    
    @test ginv ≈ [-1.0/(1.0-2/5) 0.0 0.0 0.0; 0.0 1.0-2/5 0.0 0.0; 0.0 0.0 1/25 0.0; 0.0 0.0 0.0 1/(25sin(π/3)^2)]
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

end

