using Skylight, Test

@testset "Random points and vectors" begin
    
    @testset "Cartesian coordinates" begin

        N = 5
        v = zeros(3, N)
        coord_system = Skylight.CartesianKind()

        Skylight.random_uniform_points_unit_sphere!(v, coord_system)

        for i in 1:N
            @test v[1,i]^2 + v[2,i]^2 + v[3,i]^2 ≈ 1.0 atol = 1e-14
        end

        Skylight.random_uniform_points_unit_hemisphere!(v, coord_system)

        for i in 1:N
            @test v[3,i] >= 0.0
            @test v[1,i]^2 + v[2,i]^2 + v[3,i]^2 ≈ 1.0 atol = 1e-14
        end

        Skylight.random_uniform_points_unit_hemisphere_xaxis!(v, coord_system)

        for i in 1:N
            @test v[1,i] >= 0.0
            @test v[1,i]^2 + v[2,i]^2 + v[3,i]^2 ≈ 1.0 atol = 1e-14
        end
        
        angular_radius_in_degrees = 45
        Skylight.random_uniform_points_unit_spherical_cap!(v, angular_radius_in_degrees, coord_system)
        
        for i in 1:N
            @test v[3,i] >= 1/sqrt(2)
            @test v[1,i]^2 + v[2,i]^2 + v[3,i]^2 ≈ 1.0 atol = 1e-14
        end        

        radius = 5.0

        Skylight.random_uniform_points_disk!(v, radius, coord_system)
        
        for i in 1:N
            @test v[1,i]^2 + v[2,i]^2 <= radius^2 
        end

        inner_radius = 5.0
        outer_radius = 10.0

        Skylight.random_uniform_points_annulus!(v, inner_radius, outer_radius, coord_system)

        for i in 1:N
            @test inner_radius^2 <= v[1,i]^2 + v[2,i]^2 <= outer_radius^2 
        end
    
    end

    @testset "Spherical coordinates" begin

        N = 5
        v = zeros(3, N)
        coord_system = Skylight.SphericalKind()

        Skylight.random_uniform_points_unit_sphere!(v, coord_system)

        for i in 1:N
            @test v[1,i] == 1.0 
        end
        
        Skylight.random_uniform_points_unit_hemisphere!(v, coord_system)

        for i in 1:N
            @test v[1,i] == 1.0 
            @test v[2,i] <= π/2
        end
        
        angular_radius_in_degrees = 45
        Skylight.random_uniform_points_unit_spherical_cap!(v, angular_radius_in_degrees, coord_system)
        
        for i in 1:N
            @test v[1,i] == 1.0 
            @test v[2,i] <= π/4
        end

        radius = 5.0
        Skylight.random_uniform_points_disk!(v, radius, coord_system)
        
        for i in 1:N
            @test v[1,i] <= radius
        end
        
        inner_radius = 5.0
        outer_radius = 10.0

        Skylight.random_uniform_points_annulus!(v, inner_radius, outer_radius, coord_system)

        for i in 1:N
            @test inner_radius <= v[1,i] <= outer_radius
        end
    
    end

end