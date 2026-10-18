using Skylight, Test, LinearAlgebra

# Independent reference roots from the corrected Kerr-(A)dS metric and
# hydrostatic enthalpy integral (M=1, prograde, lambda=.7).
const TORUS_REFERENCES = (
    (.1,0.,3.80349195826784,4.01007254465184,8.92491973976070,34.9164437270402,6.8339694335607e8),
    (.1,-.001,3.80076060985184,4.01895346738573,8.27044961405388,18.3247247172941,6.1073392808713e8),
    (.99,0.,2.17930534817025,1.23674442817480,1.87497794585612,5.81916324549282,3.29987880478673e9),
    (.99,-.001,2.17790913125547,1.23547506709927,1.87048285569041,5.71325683814677,3.29841419811861e9))

@testset "Four equilibrium tori" begin
    for (a,R0,l0,cusp,center,outer,Tec) in TORUS_REFERENCES
        s = FRKerrSpacetime(M=1.,a=a,R0=R0)
        model = IonTorus(s; λ=.7,ϵc=1e-17,n=1.5,β=.5,ξ=.01,
            radiative_process=SynchrotronAndBremsstrahlung())
        @test model.l0 ≈ l0 rtol=2e-9
        @test model.rcusp ≈ cusp rtol=2e-8
        @test model.rcenter ≈ center rtol=2e-8
        @test model.router ≈ outer rtol=2e-8
        @test model.Tec ≈ Tec rtol=2e-8
        @test model.Hc ≈ model.potential_at_center-model.potential_at_surface
        @test event_horizon_radius(s)+.001 < model.rcusp < model.rcenter < model.router < 50
        x = [0.,model.rcenter,π/2,0.]
        g = metric(x,s)
        u = zeros(4)
        rest_frame_four_velocity!(u,x,g,s,model,SphericalTopology())
        @test u[1]>0 && u[4]>0
        @test norm_squared(u,g) ≈ -1 atol=2e-12
        ucov = g*u
        @test -ucov[4]/ucov[1] ≈ model.l0 rtol=1e-12
        @test Skylight.energy_density(x,s,model) ≈ model.ϵc rtol=1e-12
        @test Skylight.electron_temperature(x,s,model) ≈ model.Tec rtol=1e-12
        @test Skylight.electron_temperature(x,s,model)/Skylight.ion_temperature(x,s,model) ≈ model.ξ
        @test Skylight.magnetic_field(x,s,model) ≈ Skylight.magnetic_field(model.ϵc,model)
        @test model.μₑ ≈ 8/7
        @test model.μᵢ ≈ 16/13
        for ω in (0.1,0.5,0.9)
            ρ = Skylight.energy_density(ω,model)
            p = Skylight.pressure(ρ,model)
            @test (model.n+1)*log1p(p/(ρ*PhysicalConstants.c^2)) ≈ ω*model.Hc rtol=1e-12
        end
        # Differential Euler equation, tested independently at an intermediate radius.
        r = (model.rcenter+model.router)/2
        h = 1e-5*r
        xp,xm = [0.,r+h,π/2,0.],[0.,r-h,π/2,0.]
        dp = (Skylight.pressure(xp,s,model)-Skylight.pressure(xm,s,model))/(2h)
        dW = (Skylight.torus_potential(xp,s,model)-Skylight.torus_potential(xm,s,model))/(2h)
        ρ = Skylight.energy_density([0.,r,π/2,0.],s,model)
        @test dp/(ρ*PhysicalConstants.c^2+Skylight.pressure(ρ,model)) ≈ dW rtol=2e-7
        for outside in ([0.,.9model.rcusp,π/2,0.],[0.,1.1model.router,π/2,0.],
            [0.,model.rcenter,0.,0.])
            @test Skylight.energy_density(outside,s,model)==0
            @test Skylight.electron_temperature(outside,s,model)==0
            j,alpha=ones(3),ones(3)
            E=PhysicalConstants.h .* [86e9,230e9,345e9]
            go=metric(outside,s)
            Skylight.rest_frame_emissivity!(j,outside,E,go,s,model,SphericalTopology(),nothing,nothing)
            Skylight.rest_frame_absorptivity!(alpha,outside,E,go,s,model,SphericalTopology(),nothing,nothing)
            @test all(iszero,j) && all(iszero,alpha)
        end
        E=PhysicalConstants.h .* [86e9,230e9,345e9]
        j,alpha=zeros(3),zeros(3)
        Skylight.rest_frame_emissivity!(j,x,E,g,s,model,SphericalTopology(),nothing,nothing)
        Skylight.rest_frame_absorptivity!(alpha,x,E,g,s,model,SphericalTopology(),nothing,nothing)
        @test all(>(0),j) && all(>(0),alpha)
        @test j ≈ alpha .* planck_function.(E,model.Tec) rtol=2e-12
    end
end

@testset "Explicit thermal prescription and torus domain" begin
    s=KerrSpacetimeBoyerLindquistCoordinates(M=1.,a=.5)
    @test_throws ArgumentError IonTorus(s; Tec=1e10,λ=.7,β=.5,ξ=.01)
    model=IonTorus(s; Tec=1e10,thermodynamics=:prescribed,λ=.7,β=.5,ξ=.01)
    @test model.Hc ≈ 0.353638462543 rtol=1e-10
    @test Skylight.energy_density(.5,model)/model.ϵc ≈ .334973352197 rtol=1e-10
    @test Skylight.energy_density(-.1,model)==0
    @test Skylight.energy_density(0.,model)==0
    @test Skylight.energy_density(1.,model) ≈ model.ϵc
    @test_throws AssertionError IonTorus(s; λ=0.)
    @test_throws AssertionError IonTorus(s; λ=1.)
    @test_throws ArgumentError IonTorus(KerrSpacetimeKerrSchildCoordinates(M=1.,a=.5))
    @test_throws ArgumentError IonTorus(FRKerrSpacetime(M=1.,a=.5,R0=.001))
    retro=IonTorus(s; λ=.7,rotation_sense=RetrogradeRotation())
    @test retro.rcenter>10
    @test retro.rcusp<retro.rcenter<retro.router
    @test retro.l0<0
    @test Skylight.charge_squared_ion_density(model.ϵc,model) ≈ model.ϵc/PhysicalConstants.mu
end
