using Skylight, Test, StaticArrays

# A frequency-independent homogeneous medium has an analytic transfer solution.
# Coefficients are in CGS; the ray parameter and positions are in code units.
struct UniformSlab <: AbstractRadiativeModel
    j::Float64
    alpha::Float64
    velocity::Float64
end

@testset "Failed rays are rejected" begin
    problem = ODEProblem((u,p,t)->u, [1.0], (0.0,1.0))
    good = solve(problem,Tsit5())
    @test Skylight.checked_ray_output(good,1) == (good.u[end],false)
    failed = solve(problem,Tsit5(); maxiters=1,dt=1e-3,adaptive=false,verbose=false)
    @test_throws ErrorException Skylight.checked_ray_output(failed,2)
end
Skylight.isvacuum(::UniformSlab) = Skylight.NonVacuum()
function Skylight.rest_frame_four_velocity!(v, x, g, s, m::UniformSlab, topology)
    v .= (1.0, m.velocity, 0.0, 0.0)
    normalize_timelike!(v, g)
end
function Skylight.rest_frame_emissivity!(j, x, e, g, s, m::UniformSlab, topology, sc, mc)
    fill!(j, m.j)
end
function Skylight.rest_frame_absorptivity!(alpha, x, e, g, s, m::UniformSlab, topology, sc, mc)
    fill!(alpha, m.alpha)
end

function slab_solution(; mass_unit=1.0, alpha=2e-6, velocity=0.0, affine_scale=1.0)
    spacetime = MinkowskiSpacetimeCartesianCoordinates()
    camera = PinholeCamera(position=[0.0, 10.0, 0.0, 0.0],
        horizontal_aperture_in_degrees=1.0, vertical_aperture_in_degrees=1.0,
        horizontal_number_of_pixels=1, vertical_number_of_pixels=1)
    model = UniformSlab(7e-30, alpha, velocity)
    config = NonVacuumOTEConfigurations(spacetime=spacetime, camera=camera,
        radiative_model=model, unit_mass_in_solar_masses=mass_unit,
        observation_energies=[1e-10, 3e-10, 3e-8] ./ affine_scale)
    initial = zeros(14, 1)
    initial[1:4, 1] .= camera.position
    initial[5:8, 1] .= affine_scale .* [-1.0, -1.0, 0.0, 0.0]
    cache = Skylight.transfer_cache(config, nothing, Inf)
    D = 2.0
    problem = ODEProblem(Skylight.non_vacuum_equations, SVector{14}(initial[:, 1]),
        (0.0, D / affine_scale), cache)
    solution = solve(problem, Tsit5(); reltol=1e-10, abstol=1e-13,
        save_everystep=false)
    output = reshape(collect(solution.u[end]), :, 1)
    intensity = observed_specific_intensities(initial, output, config)[:, 1]
    L = geometrized_to_CGS(1.0, Dimensions.length, config)
    q = (1 - velocity) / sqrt(1 - velocity^2)
    tau = alpha * L * q * D
    expected = alpha == 0 ? model.j * L * D / q^2 :
        model.j / (alpha * q^3) * (-expm1(-tau))
    return (; tau=output[9:11, 1], intensity, expected, expected_tau=tau)
end

@testset "CGS transfer along normalized geodesics" begin
    for mass_unit in (1.0, 7.0), alpha in (0.0, 2e-6), velocity in (0.0, 0.2)
        result = slab_solution(; mass_unit, alpha, velocity)
        @test result.tau ≈ fill(result.expected_tau, 3) rtol=1e-9 atol=1e-13
        @test result.intensity ≈ fill(result.expected, 3) rtol=2e-8 atol=0
    end
    original = slab_solution()
    rescaled = slab_solution(affine_scale=3.0)
    @test rescaled.tau ≈ original.tau rtol=1e-9
    @test rescaled.intensity ≈ original.intensity rtol=2e-8
    thin = slab_solution(alpha=0.0)
    larger = slab_solution(alpha=0.0, mass_unit=4.0)
    @test larger.intensity ≈ 4 .* thin.intensity rtol=2e-8
end
