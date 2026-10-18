using Skylight, Test, LinearAlgebra

@testset "FRKerr Kerr limit, dimensions and orbit branches" begin
    for M in (1.0, 2.3), astar in (0.0, 0.1, 0.5, 0.99)
        s = FRKerrSpacetime(M=M,a=M*astar,R0=0.0)
        kerr = KerrSpacetimeBoyerLindquistCoordinates(M=M,a=M*astar)
        @test event_horizon_radius(s) ≈ M*(1+sqrt(1-astar^2))
        @test event_horizon_radius(kerr) ≈ event_horizon_radius(s)
        @test last(horizons(s)) ≈ event_horizon_radius(s)
        for theta in (0.4, π/2, 2.1)
            x = [0.0,8M,theta,0.0]
            @test metric(x,s) ≈ metric(x,kerr) rtol=1e-13
            @test christoffel(x,s) ≈ christoffel(x,kerr) rtol=1e-11 atol=1e-13
        end
        for sense in (ProgradeRotation(),RetrogradeRotation())
            @test isco_radius(s,sense) ≈ isco_radius(kerr,sense) rtol=1e-9
            @test mbco_radius(s,sense) ≈ mbco_radius(kerr,sense) rtol=1e-9
            x = [0.0,12M,π/2,0.0]
            @test circular_geodesic_specific_angular_momentum(x,s,sense) ≈
                circular_geodesic_specific_angular_momentum(x,kerr,sense) rtol=1e-12
        end
    end
    base = FRKerrSpacetime(M=1.0,a=0.1,R0=-0.001)
    scaled = FRKerrSpacetime(M=3.0,a=0.3,R0=-0.001/9)
    sense = ProgradeRotation()
    @test isco_radius(scaled,sense) ≈ 3isco_radius(base,sense) rtol=1e-9
    @test mbco_radius(scaled,sense) ≈ 3mbco_radius(base,sense) rtol=1e-9
    @test event_horizon_radius(scaled) ≈ 3event_horizon_radius(base) rtol=1e-12
    two_roots = FRKerrSpacetime(M=1.,a=.1,R0=-.01)
    @test 4.02 < mbco_radius(two_roots,sense) < 4.04
    @test mbco_radius(two_roots,sense) < isco_radius(two_roots,sense)
    @test Skylight.outermost_circular_orbit_radius(base,sense) == Inf
    ds = FRKerrSpacetime(M=1.,a=.99,R0=.001)
    @test event_horizon_radius(ds) < cosmologic_horizon_radius(ds)
    @test cosmologic_horizon_radius(ds) ≈ 108.53056 rtol=1e-7
    @test_throws ArgumentError mbco_radius(ds,sense)
    @test_throws ArgumentError cosmologic_horizon_radius(base)
    @test_throws ArgumentError event_horizon_radius(FRKerrSpacetime(M=1.,a=.99,R0=.7))
    @test_throws ArgumentError event_horizon_radius(FRKerrSpacetime(M=1.,a=.99,R0=-.2))
end

@testset "Constant-curvature Einstein geometry" begin
    # Finite-difference derivatives of independently evaluated connections test
    # R_mn=(R0/4)g_mn off the equator, beyond just reproducing Kerr.
    for R0 in (-.001,.001)
        s = FRKerrSpacetime(M=1.,a=.5,R0=R0)
        x = [0.,8.,1.1,0.]
        G = christoffel(x,s)
        dG = zeros(4,4,4,4)
        for d in (2,3)
            h = d==2 ? 1e-4 : 1e-5
            plus,minus=copy(x),copy(x)
            plus[d]+=h; minus[d]-=h
            dG[:,:,:,d] .= (christoffel(plus,s)-christoffel(minus,s))/(2h)
        end
        ricci=zeros(4,4)
        for i in 1:4, j in 1:4, k in 1:4
            ricci[i,j] += dG[k,i,j,k]-dG[k,i,k,j]
            for l in 1:4
                ricci[i,j] += G[k,k,l]*G[l,i,j]-G[k,j,l]*G[l,i,k]
            end
        end
        @test ricci ≈ R0/4 .* metric(x,s) rtol=2e-6 atol=2e-8
    end
end
