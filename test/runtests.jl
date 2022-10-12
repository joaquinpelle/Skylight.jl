using Skylight, Test, SafeTestsets

@time begin

    @time @safetestset "Geometry" begin include("unit/utils/geometry_test.jl") end
    @time @safetestset "Random points" begin include("unit/utils/randompoints_test.jl") end
    @time @safetestset "Tetrads" begin include("unit/utils/tetrads_test.jl") end
    @time @safetestset "Emission models" begin include("unit/emissionmodels_test.jl") end
    @time @safetestset "Spacetimes" begin include("unit/spacetimes_test.jl") end
    @time @safetestset "Initial data cache" begin include("unit/initialdata/cache_test.jl") end
    @time @safetestset "Image plane" begin include("unit/initialdata/imageplane_test.jl") end
    @time @safetestset "Pixel position and momentum" begin include("unit/initialdata/pixel_position_momentum_test.jl") end
    @time @safetestset "Observer to emitter" begin include("unit/initialdata/observertoemitter_test.jl") end
    @time @safetestset "Emitter to observer" begin include("unit/initialdata/emittertoobserver_test.jl") end

end