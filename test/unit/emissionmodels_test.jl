using Skylight, Test

@testset "Synthetic polar cap" begin
    
    par = Skylight.SyntheticPolarCapParameters(number_of_points=10, 
                                               NS_radius=5.0, 
                                               misalignment_angle_in_degrees=90,
                                               angular_radius_in_degrees=60, 
                                               temperature=rand())

    @test par.misalignment_angle_in_radians ≈ π/2
    @test par.angular_radius_in_radians ≈ π/3
    @test Skylight.get_number_of_points(par) == 10

    dataframe = Skylight.synthetic_polar_cap(par)

    @test size(dataframe,2) == 10
    @test dataframe[4,:] == [par.temperature for i in 1:10]

end