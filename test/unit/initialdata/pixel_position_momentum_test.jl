using Skylight, Test

@testset "position" begin

    @testset "equatorial plane" begin

        coord_system = Skylight.SphericalClass()
            
        image_plane = ImagePlane(distance = 1.0,
                                observer_inclination_in_degrees = 90.0,
                                horizontal_side_image_plane = 1.0,
                                vertical_side_image_plane = 1.0,
                                horizontal_number_of_nodes = 3,
                                vertical_number_of_nodes = 3)
        
        pixel_coordinates = (1.0,1.0)

        r,θ,φ = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)

        @test r ≈ sqrt(3)
        @test θ ≈ acos(1/sqrt(3))
        @test φ ≈ π/4

        position_spherical = [r,θ,φ]

        coord_system = Skylight.CartesianClass()
        position_cartesian = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)

        @test position_cartesian[1] ≈ 1.0
        @test position_cartesian[2] ≈ 1.0
        @test position_cartesian[3] ≈ 1.0

        @test Skylight.spherical_from_cartesian(position_cartesian) ≈ position_spherical
    end

    @testset "z axis" begin

        coord_system = Skylight.SphericalClass()
            
        image_plane = ImagePlane(distance = 1.0,
                                observer_inclination_in_degrees = 0.0,
                                horizontal_side_image_plane = 1.0,
                                vertical_side_image_plane = 1.0,
                                horizontal_number_of_nodes = 3,
                                vertical_number_of_nodes = 3)
        
        pixel_coordinates = (1.0,1.0)

        r,θ,φ = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)

        @test r ≈ sqrt(3)
        @test θ ≈ acos(1/sqrt(3))
        @test φ ≈ atan(1,-1)

        position_spherical = [r,θ,φ]

        coord_system = Skylight.CartesianClass()
        position_cartesian = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)

        @test position_cartesian[1] ≈ -1.0
        @test position_cartesian[2] ≈  1.0
        @test position_cartesian[3] ≈  1.0

        @test Skylight.spherical_from_cartesian(position_cartesian) ≈ position_spherical
    end
    
end


@testset "momentum" begin
    
    @testset "z axis" begin
        
        coord_system = Skylight.CartesianClass()
            
        image_plane = ImagePlane(distance = 1.0,
                                observer_inclination_in_degrees = 0.0,
                                horizontal_side_image_plane = 1.0,
                                vertical_side_image_plane = 1.0,
                                horizontal_number_of_nodes = 3,
                                vertical_number_of_nodes = 3)
        
                                
        pixel_coordinates = (1.0,1.0)
        
        space_position = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)
        space_momentum = Skylight.get_space_momentum_from(pixel_coordinates,image_plane,coord_system)
        
        @test space_momentum[1] ≈ 0.0 atol=1e-15
        @test space_momentum[2] ≈ 0.0
        @test space_momentum[3] ≈ 1.0 atol=1e-15
        
        center_pixel_coordinates = (0.0,0.0)
        center_space_position = Skylight.get_space_position_from(center_pixel_coordinates,image_plane,coord_system)

        @test space_momentum'*(space_position-center_space_position) ≈ 0.0  atol=1e-15

        coord_system = Skylight.SphericalClass()
        
        space_position = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)
        space_momentum = Skylight.get_space_momentum_from(pixel_coordinates,image_plane,coord_system)

        @test space_momentum[1] ≈  1.0/sqrt(3)
        @test space_momentum[2] ≈ -sqrt(2)/3
        @test space_momentum[3] ≈  0.0
        
    end

    @testset "equatorial plane" begin
        coord_system = Skylight.CartesianClass()
            
        image_plane = ImagePlane(distance = 1.0,
                                observer_inclination_in_degrees = 90.0,
                                horizontal_side_image_plane = 1.0,
                                vertical_side_image_plane = 1.0,
                                horizontal_number_of_nodes = 3,
                                vertical_number_of_nodes = 3)
        
                                
        pixel_coordinates = (1.0,1.0)
        
        space_position = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)
        space_momentum = Skylight.get_space_momentum_from(pixel_coordinates,image_plane,coord_system)
        
        @test space_momentum[1] ≈  1.0
        @test space_momentum[2] ≈  0.0
        @test space_momentum[3] ≈  0.0 atol=1e-15
        
        center_pixel_coordinates = (0.0,0.0)
        center_space_position = Skylight.get_space_position_from(center_pixel_coordinates,image_plane,coord_system)

        @test space_momentum'*(space_position-center_space_position) ≈ 0.0  atol=1e-15

        coord_system = Skylight.SphericalClass()
        
        space_position = Skylight.get_space_position_from(pixel_coordinates,image_plane,coord_system)
        space_momentum = Skylight.get_space_momentum_from(pixel_coordinates,image_plane,coord_system)

        @test space_momentum[1] ≈ 1.0/sqrt(3)
        @test space_momentum[2] ≈ 1/(3*sqrt(2))
        @test space_momentum[3] ≈ -0.5
        
    end
end
