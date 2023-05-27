using Test

@testset "bin_values_and_sum_weights tests" begin
    bins = [0, 1, 2, 3]
    values = [0.5, 0.8, 1.5, 2.5]
    weights = [1, 2, 3, 4]

    @test Skylight.bin_values_and_sum_weights(bins, values, weights) == [3, 3, 4]

    bins = [-1, 0, 1, 2]
    values = [-0.5, 0.5, 1.5, 2.5]
    weights = [1, 2, 3, 4]

    @test Skylight.bin_values_and_sum_weights(bins, values, weights) == [1, 2, 3]

    # @test_throws ArgumentError Skylight.bin_values_and_sum_weights([0,1,2], [0,1], [0,1,2])
end


@testset "create_bins tests" begin
    @test all(Skylight.create_bins(bin_size=1.0, start=0.0, stop=3.0) .== [0.0, 1.0, 2.0, 3.0])
    @test all(Skylight.create_bins(num_bins=3, start=0.0, stop=3.0) .== [0.0, 1.0, 2.0, 3.0])
    @test_throws ErrorException Skylight.create_bins(start=0.0, stop=3.0)
end
