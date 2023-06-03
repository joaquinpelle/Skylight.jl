
using Skylight, Test

# Define the macros and struct here...

@testset "@nosave and @with_nosave" begin
    @testset "Regular fields" begin
        m = AccretionDiskWithTabulatedTemperature(inner_radius=1.0, outer_radius=2.0, filename="temp.txt")
        @test m.inner_radius == 1.0
        @test m.outer_radius == 2.0
    end

    @testset "@nosave fields" begin
        m = AccretionDiskWithTabulatedTemperature(inner_radius=1.0, outer_radius=2.0, filename="temp.txt")
        @test getfield(m,:temperature_interpolator) isa NoSaveField
        @test m.temperature_interpolator isa Function  # assuming build_interpolator returns a function
    end
end
