using Skylight, Test, SafeTestsets

@time begin

    @time @safetestset "Geometry" begin include("geometry_test.jl") end
    @time @safetestset "Spacetimes" begin include("spacetimes_test.jl") end
    @time @safetestset "Image plane" begin include("image_plane_configs_test.jl") end
    @time @safetestset "Pixel position and momentum" begin include("pixel_position_momentum_test.jl") end
    @time @safetestset "Initialization OTE" begin include("OTEinitialization_test.jl") end

end