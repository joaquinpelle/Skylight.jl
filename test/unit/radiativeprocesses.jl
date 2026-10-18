using Skylight, Test, SpecialFunctions

@testset "Thermal bremsstrahlung limits and Kirchhoff law" begin
    @test Skylight.bremsstrahlung_Fee(0.)==0
    @test Skylight.bremsstrahlung_Fee(.1) ≈ .636783357 rtol=1e-8
    @test Skylight.bremsstrahlung_Fee(.99) ≈ 32.8822037 rtol=1e-8
    for theta in exp10.(range(-8,3;length=100))
        @test Skylight.bremsstrahlung_Fee(theta)>0
        @test Skylight.bremsstrahlung_Fei(theta)>0
    end
    @test Skylight.bremsstrahlung_gaunt(.04)/Skylight.bremsstrahlung_gaunt(.01) ≈ 2
    @test Skylight.bremsstrahlung_gaunt(1.) ≈ sqrt(3)/π*(log(4)-Base.MathConstants.eulergamma)
    for Te in (1e5,1e8,1e10,1e12), ratio in (1e-8,.1,10.,1000.)
        E=ratio*PhysicalConstants.k_B*Te
        j=Skylight.bremsstrahlung_emissivity(E,1e6,1.2e6,Te)
        alpha=Skylight.bremsstrahlung_absorptivity(E,1e6,1.2e6,Te)
        @test isfinite(j) && j>=0 && isfinite(alpha) && alpha>=0
        @test j ≈ alpha*planck_function(E,Te) rtol=2e-12 atol=0
    end
    E=1e-10; Te=1e9
    @test Skylight.bremsstrahlung_emissivity(E,2e6,2.4e6,Te) ≈
        4Skylight.bremsstrahlung_emissivity(E,1e6,1.2e6,Te)
    @test Skylight.bremsstrahlung_emissivity(E,0.,0.,Te)==0
    @test Skylight.bremsstrahlung_emissivity(E,1e6,1e6,0.)==0
end

@testset "Synchrotron temperature domains" begin
    for process in (Synchrotron(),SynchrotronAndBremsstrahlung())
        for Te in (0.,1e5,prevfloat(5e8),5e8,1e9,3e10,3.1e10,3.2e10,1e12)
            for nu in (86e9,230e9,345e9,1e18)
                E=PhysicalConstants.h*nu
                j=Skylight.synchrotron_emissivity(E,1e6,Te,10.,process)
                @test isfinite(j) && j>=0
                if Te>=5e8
                    alpha=Skylight.thermal_absorptivity(j,E,Te)
                    @test isfinite(alpha) && alpha>=0
                    @test alpha*planck_function(E,Te) ≈ j rtol=2e-12 atol=0
                else
                    @test j==0 #Documented limit of the fit, no invalid interpolation.
                end
            end
        end
        @test Skylight.synchrotron_emissivity(1e-15,1e6,1e9,0.,process)==0
    end
    nu0=1e8; theta=2.; ne=1e6; B=30.; R=1e12; K2=besselk(2,1/theta)
    nuc=Skylight.critical_frequency(nu0,ne,theta,B,K2,R)
    @test isfinite(nuc) && nuc>0
    x=nuc/(1.5nu0*theta^2)
    @test exp(1.8899cbrt(x)) ≈ 2.49e-10*4π*ne*R/(B*theta^3*K2)*
        (x^(-7/6)+.4x^(-17/12)+.5316x^(-5/3)) rtol=1e-10
end

@testset "Planck and Rayleigh-Jeans intensities" begin
    Te=1e7; E=1e-10*PhysicalConstants.k_B*Te
    @test planck_function(E,Te) ≈ Skylight.rayleigh_jeans_function(E,Te) rtol=1e-9
    @test planck_function(0.,Te)==0
    @test planck_function(E,0.)==0
end
