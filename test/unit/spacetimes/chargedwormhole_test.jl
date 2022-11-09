@testset "Parameters" begin
        
    @test_throws AssertionError ChargedWormholeSpacetimeSphericalCoordinates(b0=-1.0, Q=0.5)
    @test_throws AssertionError ChargedWormholeSpacetimeSphericalCoordinates(b0=1.0, Q=1.5)
    @test_throws AssertionError ChargedWormholeSpacetimeRegularCoordinates(b0=-1.0, Q=0.5)
    @test_throws AssertionError ChargedWormholeSpacetimeRegularCoordinates(b0=1.0, Q=1.5)

end

@testset "Spherical coordinates" begin

    spacetime = ChargedWormholeSpacetimeSphericalCoordinates(b0=1.0, Q=0.5)
    
    @test Skylight.coordinate_system_class(spacetime) == Skylight.SphericalClass()

    point = [rand(),5.0,π/3,0.0]

    g = zeros(4,4)
    Skylight.set_metric!(g,point,spacetime)
    
    b = 0.2
    gtt = -(1+0.25/25)
    grr = 1/(1-0.2*0.2+0.25/25)    
    
    @test g ≈ [gtt 0.0 0.0 0.0; 0.0 grr 0.0 0.0; 0.0 0.0 25.0 0.0; 0.0 0.0 0.0 25sin(π/3)^2]
    
    ginv = zeros(4,4)
    Skylight.set_metric_inverse!(ginv,point,spacetime)
    
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

end

@testset "Regular coordinates" begin

    spacetime = ChargedWormholeSpacetimeRegularCoordinates(b0=1.0, Q=0.5)
    
    @test Skylight.coordinate_system_class(spacetime) == Skylight.SphericalClass()

    point = [rand(),5.0,π/3,0.0]

    g = zeros(4,4)
    Skylight.set_metric!(g,point,spacetime)
    
    gtt = -(1+0.25/25.75)
    
    @test g ≈ [gtt 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 25.75 0.0; 0.0 0.0 0.0 25.75*sin(π/3)^2]
    
    ginv = zeros(4,4)
    Skylight.set_metric_inverse!(ginv,point,spacetime)
    
    @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

end