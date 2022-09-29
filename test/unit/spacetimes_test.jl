using Skylight, Test

@testset "metrics" begin

    @testset "minkowski metric cartesian coordinates" begin
        

        g = zeros(4,4)
        pars = MinkowskiSpacetimeParameters()
        Skylight.minkowski_metric_cartesian_coordinates!(g,rand(4),pars)
        @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

        container = zeros(4,5)
        @views g1 = container[:,1:4]
        Skylight.minkowski_metric_cartesian_coordinates!(g1,rand(4),pars)
        @test g1 == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    
    end

    @testset "minkowski metric spherical coordinates" begin
        g = zeros(4,4)
        pars = MinkowskiSpacetimeParameters()
        position = [rand(),2.0,π/2,rand()]
        Skylight.minkowski_metric_spherical_coordinates!(g,position,pars)
        @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 4.0 0.0; 0.0 0.0 0.0 4.0]

    end

    @testset "kerr metric" begin
        
        spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,
                                                                                            a=0.0))
        metric! = spacetime.metric!

        point = [rand(),1.0,0.0,0.0]

        g1 = zeros(4,4)
        metric!(g1,point,spacetime.parameters)

        @test g1 == [1.0 2.0 0.0 0.0; 2.0 3.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

        spacetime2 = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,a=1.0))
        point = [rand(),1.0,1.0,1.0]

        g2 = zeros(4,4)
        metric!(g2,point,spacetime2.parameters)

        r2 = 1.0 + sqrt(2.0)
        r = sqrt(r2)
        H2 = 2.0*r/(r2 + 1.0/r2)
        
        l = zeros(4)
        l[1] = 1.
        l[2] = (r + 1.0)/(r2 + 1.0)
        l[3] = (r - 1.0)/(r2 + 1.0)
        l[4] = 1.0/r
        
        g = zeros(4,4)
        g[1,1]=-1. + H2 * l[1]*l[1]
        g[1,2]= 0. + H2 * l[1]*l[2]
        g[1,3]= 0. + H2 * l[1]*l[3]
        g[1,4]= 0. + H2 * l[1]*l[4]
        g[2,1]= g[1,2]
        g[2,2]= 1. + H2 * l[2]*l[2]
        g[2,3]= 0. + H2 * l[2]*l[3]
        g[2,4]= 0. + H2 * l[2]*l[4]
        g[3,1]= g[1,3]
        g[3,2]= g[2,3]
        g[3,3]= 1. + H2 * l[3]*l[3]
        g[3,4]= 0. + H2 * l[3]*l[4]
        g[4,1]= g[1,4]
        g[4,2]= g[2,4]
        g[4,3]= g[3,4]
        g[4,4]= 1. + H2 * l[4]*l[4]

        @test g2 == g
        
    end
end

@testset "spacetime parameters" begin
    
    @testset "kerr parameters" begin
        
        @test_throws AssertionError KerrSpacetimeParameters(M=-1.0,a=0.0)
        @test_throws AssertionError KerrSpacetimeParameters(M=1.0,a=1.5)
    
    end
end