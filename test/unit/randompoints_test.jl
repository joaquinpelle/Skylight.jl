using Skylight, Test

@testset "Random points and vectors" begin
    
    @testset "Cartesian coordinates" begin

        N = 100
        v = zeros(3, N)
        coord_system = Skylight.CartesianKind()

        Skylight.random_uniform_points_unit_sphere!(v, coord_system)

        @test sum(v[1,:].^2 + v[2,:].^2 + v[3,:].^2 .- 1.0) ≈ 0.0  atol = 1e-14

        Skylight.random_uniform_points_unit_hemisphere!(v, coord_system)

        @test sum(v[1,:] .>= 0.0) == N 
        @test sum(v[1,:].^2 + v[2,:].^2 + v[3,:].^2 .- 1.0) ≈ 0.0  atol = 1e-14

        
        angular_radius_in_degrees = 45
        Skylight.random_uniform_points_unit_spherical_cap!(v, angular_radius_in_degrees, coord_system)
        
        @test sum(v[1,:] .>= 1/sqrt(2)) == N 
        @test sum(v[1,:].^2 + v[2,:].^2 + v[3,:].^2 .- 1.0) ≈ 0.0  atol = 1e-14
        
        radius = 5.0
        Skylight.random_uniform_points_disk!(v, radius, coord_system)
        
        @test sum(v[1,:].^2 + v[2,:].^2 .<= radius^2)  == N 
        
        inner_radius = 5.0
        outer_radius = 10.0

        Skylight.random_uniform_points_annulus!(v, inner_radius, outer_radius, coord_system)

        @test sum(inner_radius .<= v[1,:].^2 + v[2,:].^2 .<= outer_radius^2)  == N 
    
    end

    @testset "Spherical coordinates" begin

        N = 100
        v = zeros(3, N)
        coord_system = Skylight.SphericalKind()

        Skylight.random_uniform_points_unit_sphere!(v, coord_system)

        @test sum(v[1,:] .== 1.0) == N 

        Skylight.random_uniform_points_unit_hemisphere!(v, coord_system)

        @test sum(v[1,:] .== 1.0) == N
        @test sum(v[2,:] .<= π/2) == N 
        
        angular_radius_in_degrees = 45
        Skylight.random_uniform_points_unit_spherical_cap!(v, angular_radius_in_degrees, coord_system)
        
        @test sum(v[1,:] .== 1.0) == N 
        @test sum(v[2,:] .<= π/4) == N 

        radius = 5.0
        Skylight.random_uniform_points_disk!(v, radius, coord_system)
        
        @test sum(v[1,:] .<= radius)  == N 
        
        inner_radius = 5.0
        outer_radius = 10.0

        Skylight.random_uniform_points_annulus!(v, inner_radius, outer_radius, coord_system)

        @test sum(inner_radius .<= v[1,:] .<= outer_radius)  == N 
    
    end

end