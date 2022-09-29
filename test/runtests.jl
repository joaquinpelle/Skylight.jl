using Skylight, Test

@testset "Skylight.jl" begin

    @testset "geometry" begin

        @testset "vector contraction" begin
    
            e1 =  [1,0,0,0]
            e2 =  [0,1,0,0]
            e3 =  [0,0,1,0]
            e4 =  [0,0,0,1]
            
            @test contract(e1,e1) == 1
            @test contract(e1,e2) == 0
            @test contract(e1,e3) == 0
            @test contract(e1,e4) == 0
            @test contract(e2,e2) == 1
            @test contract(e2,e3) == 0
            @test contract(e2,e4) == 0
            @test contract(e3,e3) == 1
            @test contract(e3,e4) == 0
            @test contract(e4,e4) == 1
    
            rand1 = rand(4)
            rand2 = rand(4)
    
            @test contract(rand1,rand2) == contract(rand2,rand1)
    
        end
    
        @testset "index lowering" begin
    
            metric = [-1 2 0 0 ; 2 1 0 0; 0 0 1 0; 0 0 0 1]
            vec1 = [1,0,0,0]
            vec2 = [0,1,0,0]
            @test lower_index(vec1,metric) == [-1, 2, 0, 0]
            @test lower_index(vec2,metric) == [ 2, 1, 0, 0]
        
        end
    
        @testset "scalar product" begin
    
            e1 =  [1,0,0,0]
            e2 =  [0,1,0,0]
            e3 =  [0,0,1,0]
            e4 =  [0,0,0,1]
    
            metric = [-1 0 0 0 ; 0 1 0 0; 0 0 1 0; 0 0 0 1]
            
            @test scalar_product(e1,e1,metric) == -1
            @test scalar_product(e1,e2,metric) == 0
            @test scalar_product(e1,e3,metric) == 0
            @test scalar_product(e1,e4,metric) == 0
            @test scalar_product(e2,e2,metric) == 1
            @test scalar_product(e2,e3,metric) == 0
            @test scalar_product(e2,e4,metric) == 0
            @test scalar_product(e3,e3,metric) == 1
            @test scalar_product(e3,e4,metric) == 0
            @test scalar_product(e4,e4,metric) == 1
    
            rand1 = rand(4)
            rand2 = rand(4)
    
            @test scalar_product(rand1,rand2,metric) == scalar_product(rand2,rand1,metric)
        
        end
    
        @testset "norm squared" begin
    
            metric = [-1 2 0 0 ; 2 1 0 0; 0 0 1 0; 0 0 0 1]
            vec1 = [1,0,0,0]
            vec2 = [0,1,0,0]
            vec3 = [0,0,1,0]
            @test norm_squared(vec1,metric) == -1
            @test norm_squared(vec2,metric) == 1
            @test norm_squared(vec3,metric) == 1
    
        end
    
        @testset "cartesian to spherical" begin
    
            vec1 = [1,0,0]
            vec2 = [0,1,0]
            vec3 = [0,1,1]
    
            @test cartesian_to_spherical(vec1) ≈ [1,π/2,0]
            @test cartesian_to_spherical(vec2) ≈ [1,π/2,π/2]    
            @test cartesian_to_spherical(vec3) ≈ [sqrt(2),π/4,π/2]
    
        end
    end
    
    @testset "metrics" begin
        @testset "minkowski metric cartesian coordinates" begin
            
            g = zeros(4,4)
            pars = MinkowskiSpacetimeParameters()
            minkowski_metric_cartesian_coordinates!(g,rand(4),pars)
            @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    
            container = zeros(4,5)
            @views g1 = container[:,1:4]
            minkowski_metric_cartesian_coordinates!(g1,rand(4),pars)
            @test g1 == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
        
        end
    
        @testset "minkowski metric spherical coordinates" begin
            g = zeros(4,4)
            pars = MinkowskiSpacetimeParameters()
            position = [rand(),2.0,π/2,rand()]
            minkowski_metric_spherical_coordinates!(g,position,pars)
            @test g == [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 4.0 0.0; 0.0 0.0 0.0 4.0]
    
        end
    
        @testset "kerr metric" begin
            
            spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,
                                                                                                a=0.0))
            metric! = spacetime.metric!
    
            point = [rand(),1.0,0.0,0.0]
    
            g1 = zeros(4,4)
            metric!(g1,point,spacetime.parameters)
    
            @test g1 == [1.0 2.0 0.0 0.0; 2.0 3.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    
            spacetime2 = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,a=1.0))
            point = [rand(),1.0,1.0,1.0]
    
            g2 = zeros(4,4)
            metric!(g2,point,spacetime2.parameters)
    
            r2 = 1.0 + sqrt(2.0)
            r = sqrt(r2)
            H2 = 2.0*r/(r2 + 1.0/r2)
            
            l = zeros(4)
            l[1] = 1.
            l[2] = (r + 1.0)/(r2 + 1.0)
            l[3] = (r - 1.0)/(r2 + 1.0)
            l[4] = 1.0/r
            
            g = zeros(4,4)
            g[1,1]=-1. + H2 * l[1]*l[1]
            g[1,2]= 0. + H2 * l[1]*l[2]
            g[1,3]= 0. + H2 * l[1]*l[3]
            g[1,4]= 0. + H2 * l[1]*l[4]
            g[2,1]= g[1,2]
            g[2,2]= 1. + H2 * l[2]*l[2]
            g[2,3]= 0. + H2 * l[2]*l[3]
            g[2,4]= 0. + H2 * l[2]*l[4]
            g[3,1]= g[1,3]
            g[3,2]= g[2,3]
            g[3,3]= 1. + H2 * l[3]*l[3]
            g[3,4]= 0. + H2 * l[3]*l[4]
            g[4,1]= g[1,4]
            g[4,2]= g[2,4]
            g[4,3]= g[3,4]
            g[4,4]= 1. + H2 * l[4]*l[4]
    
            @test g2 == g
            
        end
    end
    
    @testset "spacetime parameters" begin
        
        @testset "kerr parameters" begin
            
            @test_throws AssertionError KerrSpacetimeParameters(M=-1.0,a=0.0)
            @test_throws AssertionError KerrSpacetimeParameters(M=1.0,a=1.5)
        
        end
    end
    
    @testset "image plane" begin
    
        @testset "image plane area" begin
    
            image_plane = ImagePlane(observer_distance = 1.0,
                                        observer_inclination_in_degrees = 137.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
    
            @test pixel_area(image_plane) == 0.25
            @test area(image_plane) == 2.25
    
        end
    
        @testset "image plane" begin
            spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,a=0.5))
            
            image_plane = ImagePlane(observer_distance = 1.0,
                                        observer_inclination_in_degrees = 90.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
            
            configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                    image_plane = image_plane,
                                                    initial_times = [0.0])
            
            @test image_plane.observer_inclination_in_radians ≈ π/2
            @test pixel_area(image_plane) == 0.25
            @test area(image_plane) == 2.25
    
            iterator1 = get_pixel_coordinates(configurations)
            iterator2 = Iterators.product([-0.5, 0.0, 0.5], [-0.5, 0.0, 0.5])
    
            @test collect(iterator1) == collect(iterator2)
    
        end
    end
    
    @testset "configurations initialization" begin
        spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0,a=0.5))
    
        image_plane = ImagePlane(observer_distance = 1.0,
                                    observer_inclination_in_degrees = 137.0,
                                    horizontal_side_image_plane = 1.0,
                                    vertical_side_image_plane = 1.0,
                                    horizontal_number_of_nodes = 3,
                                    vertical_number_of_nodes = 3)
        
        configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                   image_plane = image_plane,
                                                   initial_times = [0.0,1.0])
        
        rays = zero_rays_on_grid(configurations)
        @test sum(rays) == 0.0
        @test length(rays)/8 == 18
        @test get_initial_times(configurations) == [0.0, 1.0]
    
    end
    
    @testset "space position and momentum from pixel coordinates" begin
    
        @testset "position" begin
    
            @testset "equatorial plane" begin
                coord_system = SphericalCoordinates()
                    
                image_plane = ImagePlane(observer_distance = 1.0,
                                        observer_inclination_in_degrees = 90.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
                
                pixel_coordinates = (1.0,1.0)
    
                r,θ,φ = get_space_position_from(pixel_coordinates,image_plane,coord_system)
    
                @test r ≈ sqrt(3)
                @test θ ≈ acos(1/sqrt(3))
                @test φ ≈ π/4
    
                position_spherical = [r,θ,φ]
    
                coord_system = CartesianCoordinates()
                position_cartesian = get_space_position_from(pixel_coordinates,image_plane,coord_system)
    
                @test position_cartesian[1] ≈ 1.0
                @test position_cartesian[2] ≈ 1.0
                @test position_cartesian[3] ≈ 1.0
    
                @test cartesian_to_spherical(position_cartesian) ≈ position_spherical
            end
    
            @testset "z axis" begin
    
                coord_system = SphericalCoordinates()
                    
                image_plane = ImagePlane(observer_distance = 1.0,
                                        observer_inclination_in_degrees = 0.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
                
                pixel_coordinates = (1.0,1.0)
    
                r,θ,φ = get_space_position_from(pixel_coordinates,image_plane,coord_system)
    
                @test r ≈ sqrt(3)
                @test θ ≈ acos(1/sqrt(3))
                @test φ ≈ atan(1,-1)
    
                position_spherical = [r,θ,φ]
    
                coord_system = CartesianCoordinates()
                position_cartesian = get_space_position_from(pixel_coordinates,image_plane,coord_system)
    
                @test position_cartesian[1] ≈ -1.0
                @test position_cartesian[2] ≈  1.0
                @test position_cartesian[3] ≈  1.0
    
                @test cartesian_to_spherical(position_cartesian) ≈ position_spherical
            end
            
        end
    
    
        @testset "momentum" begin
            
            @testset "z axis" begin
                
                coord_system = CartesianCoordinates()
                    
                image_plane = ImagePlane(observer_distance = 1.0,
                                        observer_inclination_in_degrees = 0.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
                
                                        
                pixel_coordinates = (1.0,1.0)
                
                space_position = get_space_position_from(pixel_coordinates,image_plane,coord_system)
                space_momentum = get_space_momentum_from(pixel_coordinates,image_plane,coord_system)
                
                @test space_momentum[1] ≈ 0.0 atol=1e-15
                @test space_momentum[2] ≈ 0.0
                @test space_momentum[3] ≈ 1.0 atol=1e-15
                
                center_pixel_coordinates = (0.0,0.0)
                center_space_position = get_space_position_from(center_pixel_coordinates,image_plane,coord_system)
    
                @test space_momentum'*(space_position-center_space_position) ≈ 0.0  atol=1e-15
    
                coord_system = SphericalCoordinates()
                
                space_position = get_space_position_from(pixel_coordinates,image_plane,coord_system)
                space_momentum = get_space_momentum_from(pixel_coordinates,image_plane,coord_system)
    
                @test space_momentum[1] ≈  1.0/sqrt(3)
                @test space_momentum[2] ≈ -sqrt(2)/3
                @test space_momentum[3] ≈  0.0
                
            end
    
            @testset "equatorial plane" begin
                coord_system = CartesianCoordinates()
                    
                image_plane = ImagePlane(observer_distance = 1.0,
                                        observer_inclination_in_degrees = 90.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
                
                                        
                pixel_coordinates = (1.0,1.0)
                
                space_position = get_space_position_from(pixel_coordinates,image_plane,coord_system)
                space_momentum = get_space_momentum_from(pixel_coordinates,image_plane,coord_system)
                
                @test space_momentum[1] ≈  1.0
                @test space_momentum[2] ≈  0.0
                @test space_momentum[3] ≈  0.0 atol=1e-15
                
                center_pixel_coordinates = (0.0,0.0)
                center_space_position = get_space_position_from(center_pixel_coordinates,image_plane,coord_system)
    
                @test space_momentum'*(space_position-center_space_position) ≈ 0.0  atol=1e-15
    
                coord_system = SphericalCoordinates()
                
                space_position = get_space_position_from(pixel_coordinates,image_plane,coord_system)
                space_momentum = get_space_momentum_from(pixel_coordinates,image_plane,coord_system)
    
                @test space_momentum[1] ≈ 1.0/sqrt(3)
                @test space_momentum[2] ≈ 1/(3*sqrt(2))
                @test space_momentum[3] ≈ -0.5
                
            end
        end
    end
    
    @testset "initialization" begin
    
        @testset "four-momentum" begin
        
            @testset "dumps in container" begin
        
                container = zeros(4,5)
    
                dump_∂t_in!(container)
                @test container == [0 0 0 0 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
    
                spacetime = MinkowskiSpacetimeCartesianCoordinates()
    
                image_plane = ImagePlane(observer_distance = 1.0,
                                            observer_inclination_in_degrees = 137.0,
                                            horizontal_side_image_plane = 1.0,
                                            vertical_side_image_plane = 1.0,
                                            horizontal_number_of_nodes = 3,
                                            vertical_number_of_nodes = 3)
                
                configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                        image_plane = image_plane,
                                                        initial_times = [0.0,1.0])
    
                position = rand(4)
                dump_metric_in!(container,position,spacetime)
                
                @test container == [-1.0 0.0 0.0 0.0 1.0; 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0]
    
            end
    
            @testset "set null ingoing past directed four-momentum" begin
    
                spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0, a=0.99))
    
                position = [rand(),5.0,rand(),rand()]
                momentum = rand(4)
                momentum[1] = 0.0
    
                container = zeros(4,5)
                dump_∂t_in!(container)
                dump_metric_in!(container,position,spacetime)
    
                set_null!(momentum,container)
                metric! = spacetime.metric!
                
                g = zeros(4,4)
                g = metric!(g,position,spacetime.parameters)
                
                @test norm_squared(momentum,g) ≈ 0.0 atol=1e-15
                
                @test momentum[1] == 1.0
                set_ingoing_past_directed!(momentum)
                @test momentum[1] == -1.0
                
                momentum = rand(4)
                momentum[1] = 0.0
    
                set_null_ingoing_past_directed!(momentum,container)
                
                @test norm_squared(momentum,g) ≈ 0.0 atol=1e-15
                @test momentum[1] == -1.0
            end
        end
    
        @testset "single ray" begin
            
            spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0, a=0.0))
    
            image_plane = ImagePlane(observer_distance = sqrt(7.0),
                                        observer_inclination_in_degrees = 90.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
            
            configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                    image_plane = image_plane,
                                                    initial_times = [0.0,rand()])
    
            initial_time = configurations.initial_times[2]
            image_plane = configurations.image_plane
            coord_system = configurations.spacetime.asymptotic_coordinate_system
    
            ray = zeros(8)
            pixel_coordinates = (1.0, 1.0)
            container = zeros(4,5)
            dump_∂t_in!(container)
    
            initialize_single!(ray,initial_time,pixel_coordinates,configurations,container)
            
            @test ray[1] ==  initial_time
            @test ray[2] ≈   sqrt(7.0)
            @test ray[3] ≈   1.0
            @test ray[4] ≈   1.0
            @test ray[5] == -1.0
            @test ray[6] ≈  -3*(-2*sqrt(7)+sqrt(69))/41
            @test ray[7] ≈   0.0
            @test ray[8] ≈   0.0  atol=1e-15
            
        end
    
        @testset "all rays" begin
            
            spacetime = KerrSpacetimeKerrSchildCoordinates(parameters = KerrSpacetimeParameters(M=1.0, a=0.0))
    
            image_plane = ImagePlane(observer_distance = 3,
                                        observer_inclination_in_degrees = 90.0,
                                        horizontal_side_image_plane = 1.0,
                                        vertical_side_image_plane = 1.0,
                                        horizontal_number_of_nodes = 3,
                                        vertical_number_of_nodes = 3)
            
            configurations = OTEInitialDataConfigurations(spacetime=spacetime,
                                                    image_plane = image_plane,
                                                    initial_times = [0.1,1.5])
    
            rays = initialize_OTE(configurations)
    
            @views ray = rays[14,:]
    
            @test ray[1] ==  1.5
            @test ray[2] ≈   3.0
            @test ray[3] ≈   0.0
            @test ray[4] ≈   0.0 atol=1e-15
            @test ray[5] == -1.0
            @test ray[6] ≈  -0.2
            @test ray[7] ≈   0.0
            @test ray[8] ≈   0.0  atol=1e-15
    
            @views ray = rays[1,:]
    
            @test ray[1] ==  0.1
            @test ray[2] ≈   3.0
            @test ray[3] ≈  -0.5
            @test ray[4] ≈  -0.5 atol=1e-15
            @test ray[5] == -1.0
            @test ray[7] ≈   0.0
            @test ray[8] ≈   0.0  atol=1e-15
            
        end
    end
end

