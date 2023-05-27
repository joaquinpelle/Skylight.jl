@testset "Parameters" begin
        
    @test_throws AssertionError KerrSpacetimeKerrSchildCoordinates(M=-1.0,a=0.0)
    @test_throws AssertionError KerrSpacetimeKerrSchildCoordinates(M=1.0,a=1.5)
    @test_throws AssertionError KerrSpacetimeBoyerLindquistCoordinates(M=-1.0,a=0.0)
    @test_throws AssertionError KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=1.5)
end

@testset "Christoffel cache" begin
    
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)

    @test typeof(allocate_christoffel_cache(spacetime)) == Skylight.KerrChristoffelCache

end

@testset "Metric" begin

    @testset "Kerr-Schild coordinates" begin
            
        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.0)

        @test coordinates_topology(spacetime) == CartesianTopology()

        point = [rand(),1.0,0.0,0.0]

        g1 = zeros(4,4)
        set_metric!(g1,point,spacetime)

        @test g1 == [1.0 2.0 0.0 0.0; 2.0 3.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

        spacetime2 = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=1.0)
        point = [rand(),1.0,1.0,1.0]

        g2 = zeros(4,4)
        set_metric!(g2,point,spacetime2)

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

        spacetime3 = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)
        position = [rand(), 4.0, 5.0, 1.0]
        
        g3 = zeros(4,4)
        set_metric!(g3,point,spacetime3)

        ginv = zeros(4,4)
        set_metric_inverse!(ginv,point,spacetime3)
        @test g3*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    end

    @testset "Boyer-Lindquist coordinates" begin
        
        spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.5)

        @test coordinates_topology(spacetime) == SphericalTopology()

        point = [rand(),5.0,π/3,0.0]

        g = zeros(4,4)
        set_metric!(g,point,spacetime)
        
        Σ = 25+0.25*0.25
        Δ = 25-10+0.25

        gtt = -1.0 + 10/Σ
        gtφ = -15/(4*Σ)
        gφφ = 0.75*(25+0.25+2.5*0.75/Σ)
        

        @test g ≈ [gtt 0.0 0.0 gtφ; 0.0 Σ/Δ 0.0 0.0; 0.0 0.0 Σ 0.0; gtφ 0.0 0.0 gφφ]
        
        ginv = zeros(4,4)
        set_metric_inverse!(ginv,point,spacetime)
        
        det = gtt*gφφ-gtφ^2

        @test ginv ≈ [gφφ/det 0.0 0.0 -gtφ/det; 0.0 Δ/Σ 0.0 0.0; 0.0 0.0 1/Σ 0.0; -gtφ/det 0.0 0.0 gtt/det]
        @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    end
end

@testset "Christoffel" begin
    
    spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.9)
    position = rand(4)
    cache = allocate_christoffel_cache(spacetime)

    Γ = zeros(4,4,4)

    set_christoffel!(Γ,position,spacetime,cache)

    @test sum([Γ[i] != 0.0 for i in 1:64]) == 64

end