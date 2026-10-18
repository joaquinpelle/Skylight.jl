using Skylight, Test, HDF5, DelimitedFiles, LinearAlgebra
include(joinpath(@__DIR__,"..","..","run","iontorus.jl"))

@testset "Four-case example and saved observables" begin
    mktempdir() do output
        for case in IonTorusExample.CASES
            result=IonTorusExample.run_case(case; pixels=16,output)
            @test size(result.intensity_nu)==(3,16*12)
            @test all(isfinite,result.final)
            @test all(>=(0),result.intensity_nu)
            @test maximum(result.intensity_nu)>0
            @test all(>=(0),result.final[9:11,:])
            # Stationary/axisymmetric null rays conserve their two Killing
            # momenta. This also checks the actual multi-threaded AD path.
            s=result.config.spacetime
            g0=metric(result.config.camera.position,s)
            initial_cov=g0*result.initial[5:8,:]
            for i in axes(result.final,2)
                k=result.final[5:8,i]
                gf=metric(result.final[1:4,i],s)
                kcov=gf*k
                @test kcov[1] ≈ initial_cov[1,i] rtol=3e-6
                @test kcov[4] ≈ initial_cov[4,i] rtol=3e-6 atol=1e-9
                @test abs(dot(k,kcov))/sum(abs.(k.*kcov)) < 3e-6
            end
            path=joinpath(output,case.id)
            h5open(path*".h5","r") do file
                @test read(file["I_nu"]) ≈ reshape(result.intensity_nu,3,16,12)
                @test read(file["I_nu"]) ≈ PhysicalConstants.h .* read(file["I_E"])
                @test read(file["F_nu"]) ≈ result.flux_nu
                @test read(attributes(file)["case_id"])==case.id
            end
            table=readdlm(path*"_spectrum.csv",','; skipstart=1)
            @test table[:,2] ≈ result.flux_nu
            @test table[:,3] ≈ IonTorusExample.FREQUENCIES_HZ .* result.flux_nu
            @test isfile(path*"_230GHz.svg")
            @test_throws ArgumentError IonTorusExample.run_case(case; pixels=16,output)
        end
        @test_throws ArgumentError IonTorusExample.main(["--case","unknown"])
        @test_throws ArgumentError IonTorusExample.configuration(first(IonTorusExample.CASES); pixels=7)
        c=IonTorusExample.configuration(last(IonTorusExample.CASES); pixels=4)
        @test_throws ArgumentError IonTorusExample.validate_domain(c; horizon_buffer=.2)
    end
end
