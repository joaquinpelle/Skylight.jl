using Skylight, Test

@testset "Vector contraction" begin

    contract = Skylight.contract

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

@testset "Index lowering" begin

    lower_index = Skylight.lower_index

    metric = [-1 2 0 0 ; 2 1 0 0; 0 0 1 0; 0 0 0 1]
    vec1 = [1,0,0,0]
    vec2 = [0,1,0,0]

    @test lower_index(vec1,metric) == [-1, 2, 0, 0]
    @test lower_index(vec2,metric) == [ 2, 1, 0, 0]

end

@testset "Scalar product" begin

    scalar_product = Skylight.scalar_product

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

@testset "Norm squared" begin

    norm_squared = Skylight.norm_squared

    metric = [-1 2 0 0 ; 2 1 0 0; 0 0 1 0; 0 0 0 1]
    vec1 = [1,0,0,0]
    vec2 = [0,1,0,0]
    vec3 = [0,0,1,0]
    @test norm_squared(vec1,metric) == -1
    @test norm_squared(vec2,metric) == 1
    @test norm_squared(vec3,metric) == 1

end

@testset "Normalization" begin

    e1 = [1.0, 0.0, 0.0, 0.0]
    e2 = [0.0, 1.0, 0.0, 0.0]

    g = [-rand() 0 0 0; 0 rand() 0 0; 0 0 rand() 0; 0 0 0 rand()]

    Skylight.normalize_timelike!(e1, g)
    Skylight.normalize_spacelike!(e2, g)

    @test Skylight.norm_squared(e1,g) ≈ -1.0
    @test Skylight.norm_squared(e2,g) ≈  1.0

end

@testset "Tangent vector zaxis rotation" begin
    
    vector = zeros(4)
    vector1 = zeros(4)
    position = rand(4)
    metric = [-rand() 0 0 0; 0 rand() 0 0; 0 0 rand() 0; 0 0 0 rand()]
    Ω = 0.05

    @testset "Cartesian coordinates" begin

        circular_motion_four_velocity!(vector1,position,Ω,metric,CartesianTopology())
        
        vector[1] =  1.0
        vector[2] = -Ω*position[3]
        vector[3] =  Ω*position[2]

        Skylight.normalize_timelike!(vector,metric)

        @test vector1 ≈ vector
    
    end

    @testset "Spherical coordinates" begin

        circular_motion_four_velocity!(vector1,position,Ω,metric,SphericalTopology())
        
        vector[1] = 1.0
        vector[2] = 0.0
        vector[3] = 0.0
        vector[4] = Ω

        Skylight.normalize_timelike!(vector,metric)

        @test vector1 ≈ vector
    
    end

end

@testset "Cartesian to spherical" begin

    spherical_from_cartesian = Skylight.spherical_from_cartesian
    
    vec1 = [1,0,0]
    vec2 = [0,1,0]
    vec3 = [0,1,1]

    @test spherical_from_cartesian(vec1) ≈ [1,π/2,0]
    @test spherical_from_cartesian(vec2) ≈ [1,π/2,π/2]    
    @test spherical_from_cartesian(vec3) ≈ [sqrt(2),π/4,π/2]

end

@testset "Rotation around y-axis" begin
    
    v = zeros(3,3)

    v[:,1] .= [1, 0, 0]
    v[:,2] .= [0, 1, 0]
    v[:,3] .= [0, 0, 1]

    Skylight.rotate_around_y_axis!(v,90)

    @test v[:,1] ≈ [0, 0, -1]
    @test v[:,2] ≈ [0, 1,  0]
    @test v[:,3] ≈ [1, 0,  0]

end