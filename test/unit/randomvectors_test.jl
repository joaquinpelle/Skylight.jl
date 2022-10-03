using Skylight, Test

@testset "Random points and vectors" begin
    
    N = 100
    v = zeros(N,3)

    Skylight.random_isotropic_unit_vectors_sphere!(v)

    @test sum(v[:,1].^2 + v[:,2].^2 + v[:,3].^2 .- 1.0) ≈ 0.0  atol = 1e-14

    Skylight.random_isotropic_unit_vectors_hemisphere!(v)

    @test sum(v[:,3] .>= 0.0) == N 
    @test sum(v[:,1].^2 + v[:,2].^2 + v[:,3].^2 .- 1.0) ≈ 0.0  atol = 1e-14

    
    angular_radius_in_degrees = 45
    Skylight.random_isotropic_unit_vectors_cone!(v,angular_radius_in_degrees)
    
    @test sum(v[:,3] .>= 1/sqrt(2)) == N 
    @test sum(v[:,1].^2 + v[:,2].^2 + v[:,3].^2 .- 1.0) ≈ 0.0  atol = 1e-14
    
    radius = 5.0
    Skylight.random_uniform_points_disk!(v,radius)
    
    @test sum(v[:,1].^2 + v[:,2].^2 .<= radius^2)  == N 
    
    inner_radius = 5.0
    outer_radius = 10.0

    Skylight.random_uniform_points_annulus!(v, inner_radius, outer_radius)

    @test sum(inner_radius .<= v[:,1].^2 + v[:,2].^2 .<= outer_radius^2)  == N 


end