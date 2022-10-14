using Skylight, Test, SafeTestsets

@time begin

    @safetestset "Geometry" begin include("unit/utils/geometry_test.jl") end
    @safetestset "Random points" begin include("unit/utils/randompoints_test.jl") end
    @safetestset "Tetrads" begin include("unit/utils/tetrads_test.jl") end
    @safetestset "Emission models" begin include("unit/emissionmodels_test.jl") end
    @safetestset "Spacetimes" begin include("unit/spacetimes/test.jl") end
    @safetestset "Initial data cache" begin include("unit/initialdata/cache_test.jl") end
    @safetestset "Image plane" begin include("unit/initialdata/imageplane_test.jl") end
    @safetestset "Pixel position and momentum" begin include("unit/initialdata/pixel_position_momentum_test.jl") end
    @safetestset "Observer to emitter" begin include("unit/initialdata/observertoemitter_test.jl") end
    @safetestset "Emitter to observer" begin include("unit/initialdata/emittertoobserver_test.jl") end
    @safetestset "Geodesics cache" begin include("unit/geodesics/cache_test.jl") end

end