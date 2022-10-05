using Skylight, Test, SafeTestsets

@time begin

    @time @safetestset "Geometry" begin include("unit/geometry_test.jl") end
    @time @safetestset "Random vectors" begin include("unit/randomvectors_test.jl") end
    @time @safetestset "Dump containers" begin include("unit/dumpcontainers_test.jl") end
    @time @safetestset "Spacetimes" begin include("unit/spacetimes_test.jl") end
    @time @safetestset "Emission models" begin include("unit/emissionmodels_test.jl") end
    @time @safetestset "Image plane" begin include("unit/image_plane_test.jl") end
    @time @safetestset "Pixel position and momentum" begin include("unit/pixel_position_momentum_test.jl") end
    @time @safetestset "Initialization OTE" begin include("unit/OTEinitialization_test.jl") end

end