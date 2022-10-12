@testset "Parameters" begin
        
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

@testset "Boyer-Lindquist coordinates" begin
        
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