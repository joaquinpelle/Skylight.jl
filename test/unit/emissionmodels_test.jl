using Skylight, Test

@testset "Synthetic polar cap" begin
    
    model = Skylight.SyntheticPolarCap(number_of_points=10, 
                                    NS_radius=5.0,
                                    angular_speed = 0.05, 
                                    misalignment_angle_in_degrees=90,
                                    angular_radius_in_degrees=60, 
                                    temperature=rand())

    @test model.misalignment_angle_in_radians ≈ π/2
    @test model.angular_radius_in_radians ≈ π/3
    @test Skylight.get_number_of_points(model) == 10

    coord_system = Skylight.CartesianKind()
    dataframe = Skylight.synthetic_polar_cap_dataframe(model, coord_system)

    @test size(dataframe,2) == 10
    @test dataframe[4,:] == [model.temperature for i in 1:10]

    vector = zeros(4)
    position = [rand(), 3.0, 0.0, 4.0]
    gμν = [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

    Skylight.set_local_four_velocity!(vector, position, gμν, model, coord_system)

    @test vector ≈ [1.0/sqrt(0.9775), 0.0, 0.15/sqrt(0.9775), 0.0]

    @test Skylight.surface_function(position, model, coord_system) ≈ 0.0
    
    df = zeros(4)
    Skylight.set_surface_differential!(df, position, model, coord_system)
    @test df == [0.0, 2*position[2], 2*position[3], 2*position[4]]


    spacetime = Skylight.KerrSpacetimeKerrSchildCoordinates(M=1.0,a=0.5)

    gcache = zeros(4,4)
    Skylight.set_metric!(gcache, position, spacetime)

    normal = zeros(4)
    Skylight.set_unit_surface_normal!(normal, position, gcache, spacetime, model, coord_system)

    tangent_vector = [0.0, -position[3], position[2], 0.0]

    Skylight.set_metric!(gcache, position, spacetime)
    @test Skylight.norm_squared(normal, gcache) ≈ 1.0
    @test Skylight.scalar_product(normal, tangent_vector, gcache) ≈ 0.0 atol=1e-16
    
    Skylight.lower_index!(normal, gcache)
    @test normal[2]/df[2] ≈ normal[4]/df[4]

end