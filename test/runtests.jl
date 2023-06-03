using Skylight, Test, SafeTestsets

@time begin

    @safetestset "Geometry" begin include("unit/utils/bins.jl") end
    @safetestset "Geometry" begin include("unit/utils/geometry.jl") end
    @safetestset "Random points" begin include("unit/utils/randompoints.jl") end
    @safetestset "Tetrads" begin include("unit/utils/tetrads.jl") end
    @safetestset "HDF5" begin include("unit/utils/save_hdf5.jl") end
    @safetestset "Emission models" begin include("unit/radiativemodels.jl") end
    @safetestset "Spacetimes" begin include("unit/spacetimes/spacetimes.jl") end
    @safetestset "Initial data cache" begin include("unit/initialdata/cache.jl") end
    @safetestset "Image plane" begin include("unit/initialdata/imageplane.jl") end
    @safetestset "Pixel position and momentum" begin include("unit/initialdata/pixel_position_momentum.jl") end
    @safetestset "Observer to emitter" begin include("unit/initialdata/observertoemitter.jl") end
    @safetestset "Emitter to observer" begin include("unit/initialdata/emittertoobserver.jl") end
    @safetestset "Geodesics cache" begin include("unit/transfer/cache.jl") end
    @safetestset "Callbacks" begin include("unit/transfer/callbacks.jl") end
    @safetestset "Problem" begin include("unit/transfer/problem.jl") end

end