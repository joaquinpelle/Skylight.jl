using Skylight, Test

@testset "Metrics" begin

    @testset "Minkowski metric cartesian coordinates" begin
        
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

    @testset "Minkowski metric spherical coordinates" begin
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

    @testset "Schwarzschild metric Kerr-Schild coordinates" begin
        
        spacetime = SchwarzschildSpacetimeKerrSchildCoordinates(M=1.0)
        @test Skylight.coordinate_system_kind(spacetime) == Skylight.CartesianKind()


        point = [rand(),0.0,3.0,4.0]

        g = zeros(4,4)
        Skylight.set_metric!(g,point,spacetime)
        
        @test g ≈ [-1.0+2/5 0.0 6/25 8/25; 0.0 1.0 0.0 0.0; 6/25 0.0 1.0+18/125 24/125; 8/25 0.0 24/125 1.0+32/125]
        
        ginv = zeros(4,4)
        Skylight.set_metric_inverse!(ginv,point,spacetime)

        @test ginv ≈ [-1.0-2/5 0.0 6/25 8/25; 0.0 1.0 0.0 0.0; 6/25 0.0 1.0-18/125 -24/125; 8/25 0.0 -24/125 1.0-32/125]
        @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    end

    @testset "Schwarzschild metric spherical coordinates" begin

        spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
        @test Skylight.coordinate_system_kind(spacetime) == Skylight.SphericalKind()


        point = [rand(),5.0,π/3,0.0]

        g = zeros(4,4)
        Skylight.set_metric!(g,point,spacetime)
        
        @test g ≈ [-(1.0-2/5) 0.0 0.0 0.0; 0.0 1.0/(1.0-2/5) 0.0 0.0; 0.0 0.0 25.0 0.0; 0.0 0.0 0.0 25sin(π/3)^2]
        
        ginv = zeros(4,4)
        Skylight.set_metric_inverse!(ginv,point,spacetime)
        
        @test ginv ≈ [-1.0/(1.0-2/5) 0.0 0.0 0.0; 0.0 1.0-2/5 0.0 0.0; 0.0 0.0 1/25 0.0; 0.0 0.0 0.0 1/(25sin(π/3)^2)]
        @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    end

    @testset "Kerr metric Kerr-Schild coordinates" begin
        
        spacetime = KerrSpacetimeKerrSchildCoordinates(M=1.0, a=0.0)

        @test Skylight.coordinate_system_kind(spacetime) == Skylight.CartesianKind()

        point = [rand(),1.0,0.0,0.0]

        g1 = zeros(4,4)
        Skylight.set_metric!(g1,point,spacetime)

        @test g1 == [1.0 2.0 0.0 0.0; 2.0 3.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

        spacetime2 = KerrSpacetimeKerrSchildCoordinates(M=1.0,a=1.0)
        point = [rand(),1.0,1.0,1.0]

        g2 = zeros(4,4)
        Skylight.set_metric!(g2,point,spacetime2)

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
        Skylight.set_metric!(g3,point,spacetime3)

        ginv = zeros(4,4)
        Skylight.set_metric_inverse!(ginv,point,spacetime3)
        @test g3*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    end

    @testset "Kerr metric Boyer-Lindquist coordinates" begin
        
        spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.5)

        @test Skylight.coordinate_system_kind(spacetime) == Skylight.SphericalKind()

        point = [rand(),5.0,π/3,0.0]

        g = zeros(4,4)
        Skylight.set_metric!(g,point,spacetime)
        
        Σ = 25+0.25*0.25
        Δ = 25-10+0.25

        gtt = -1.0 + 10/Σ
        gtφ = -15/(4*Σ)
        gφφ = 0.75*(25+0.25+2.5*0.75/Σ)
        

        @test g ≈ [gtt 0.0 0.0 gtφ; 0.0 Σ/Δ 0.0 0.0; 0.0 0.0 Σ 0.0; gtφ 0.0 0.0 gφφ]
        
        ginv = zeros(4,4)
        Skylight.set_metric_inverse!(ginv,point,spacetime)
        
        det = gtt*gφφ-gtφ^2

        @test ginv ≈ [gφφ/det 0.0 0.0 -gtφ/det; 0.0 Δ/Σ 0.0 0.0; 0.0 0.0 1/Σ 0.0; -gtφ/det 0.0 0.0 gtt/det]
        @test g*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    end

    @testset "Johannsen Spacetime" begin
        
        spacetime = JohannsenSpacetimeBoyerLindquistCoordinates(M=1.0,
                                                                a=0.5,
                                                                α13=0.2,
                                                                α22=0.1,
                                                                α52=0.3,
                                                                ϵ3=0.5)

        @test Skylight.coordinate_system_kind(spacetime) == Skylight.SphericalKind()


        r = 5.0
        θ = π/3

        point = [rand(),r,θ,rand()]

        M = spacetime.M
        a = spacetime.a
        α13 = spacetime.α13
        α22 = spacetime.α22
        α52 = spacetime.α52
        ϵ3 = spacetime.ϵ3

        r2 = r^2
        a2 = a^2
        sinθ2 = sin(θ)^2

        Δ = r2 - 2*M*r + a2
        Σ = r2 + a2*cos(θ)^2 + ϵ3*M^3/r

        A1 = 1 + α13*(M/r)^3
        A2 = 1 + α22*(M/r)^2
        A5 = 1 + α52*(M/r)^2

        C = ((r2+a2)*A1-a2*A2*sinθ2)^2

        g = zeros(4,4)

        g[1,1] = -Σ*(Δ-a2*A2^2*sinθ2)/C
        g[1,4] = -a*((r2+a2)*A1*A2-Δ)*Σ*sinθ2/C 

        g[2,2] = Σ/(Δ*A5)

        g[3,3] = Σ

        g[4,1] = g[1,4]
        g[4,4] = Σ*sinθ2*((r2+a2)^2*A1^2-a2*Δ*sinθ2)/C
        
        g1 = zeros(4,4)
        Skylight.set_metric!(g1,point,spacetime)

        @test g ≈ g1
        
        ginv = zeros(4,4)
        Skylight.set_metric_inverse!(ginv,point,spacetime)
        
        @test g1*ginv ≈ [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    
    end

end

@testset "Spacetime parameters" begin
    
    @testset "Schwarzschild" begin
        
        @test_throws AssertionError SchwarzschildSpacetimeKerrSchildCoordinates(M=-1.0)
        @test_throws AssertionError SchwarzschildSpacetimeSphericalCoordinates(M=-1.0)
    
    end

    @testset "Kerr" begin
        
        @test_throws AssertionError KerrSpacetimeKerrSchildCoordinates(M=-1.0,a=0.0)
        @test_throws AssertionError KerrSpacetimeKerrSchildCoordinates(M=1.0,a=1.5)
        @test_throws AssertionError KerrSpacetimeBoyerLindquistCoordinates(M=-1.0,a=0.0)
        @test_throws AssertionError KerrSpacetimeBoyerLindquistCoordinates(M=1.0,a=1.5)
    end

    @testset "Johannsen" begin
        
        @test_throws AssertionError JohannsenSpacetimeBoyerLindquistCoordinates(M=-1.0,
                                                                                a=0.5,
                                                                                α13=0.2,
                                                                                α22=0.1,
                                                                                α52=0.3,
                                                                                ϵ3=0.5)

        @test_throws DomainError JohannsenSpacetimeBoyerLindquistCoordinates(M=1.0,
                                                                             a=1.5,
                                                                             α13=0.2,
                                                                             α22=0.1,
                                                                             α52=0.3,
                                                                             ϵ3=0.5)


        @test_throws AssertionError JohannsenSpacetimeBoyerLindquistCoordinates(M=1.0,
                                                                                a=sqrt(0.75),
                                                                                α13=-1.6^3,
                                                                                α22=0.1,
                                                                                α52=0.3,
                                                                                ϵ3=0.5)
        @test_throws AssertionError JohannsenSpacetimeBoyerLindquistCoordinates(M=1.0,
                                                                                a=sqrt(0.75),
                                                                                α13=0.6,
                                                                                α22=-2.1^3,
                                                                                α52=0.3,
                                                                                ϵ3=0.5)
        @test_throws AssertionError JohannsenSpacetimeBoyerLindquistCoordinates(M=1.0,
                                                                                a=sqrt(0.75),
                                                                                α13=0.6,
                                                                                α22=0.1,
                                                                                α52=-2.3^3,
                                                                                ϵ3=0.5)
                                                                                
        @test_throws AssertionError JohannsenSpacetimeBoyerLindquistCoordinates(M=1.0,
                                                                                a=sqrt(0.75),
                                                                                α13=-0.6,
                                                                                α22=0.1,
                                                                                α52=0.3,
                                                                                ϵ3=-2.5^3)
                                                                                
                                                                                
                                                                                
    end

end