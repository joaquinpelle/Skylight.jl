using Skylight, HDF5, Test

# Custom type for testing
struct CustomType
    field1::Int
    field2::Float64
    field3::String
end

# Define a custom type for testing
struct MyType2
    a::Int
    b::Float64
    c::Complex{Float64}
    d::String
    e::Function
end

function to_dict(c::CustomType)
    return Dict(:field1 => c.field1, :field2 => c.field2, :field3 => c.field3)
end

spacetime = KerrSpacetimeKerrSchildCoordinates(M = 1.0, a = 0.9)

camera = ImagePlane(distance = 500.0,
    observer_inclination_in_degrees = 45,
    horizontal_side = 20.0,
    vertical_side = 20.0,
    horizontal_number_of_pixels = 3,
    vertical_number_of_pixels = 3)

model = DummyDisk(inner_radius = 6.0, outer_radius = 18.0)

configurations = VacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0)

callback, callback_parameters = callback_setup(configurations; rhorizon_bound = 2e-1)
initial_data = rand(10, 10)
output_data = rand(10, 10)

@testset "Save and load tests" begin
    mktempdir() do dir
        # 'dir' is a string containing the path to a new temporary directory
        # This directory is automatically deleted when this block exits

        # Test case 1: save_to_hdf5
        @testset "save_to_hdf5" begin
            filename = joinpath(dir, "test.h5")
            kwargs = Dict(:key1 => "value1", :key2 => "value2")
            runs = [Skylight.Run(output_data, callback, callback_parameters, kwargs)]
            save_to_hdf5(filename, configurations, initial_data, runs)
            @test isfile(filename)
            h5open(filename, "r") do file
                @test keys(file) == ["configs", "initial_data", "run_1"]
                @test length(keys(file["run_1"])) == 4
            end
        end

        # Test case 2: append_runs_to_hdf5
        @testset "append_runs_to_hdf5" begin
            filename = joinpath(dir, "test.h5")
            kwargs = Dict(:key3 => "value3", :key4 => "value4")
            new_runs = [Skylight.Run(output_data, callback, callback_parameters, kwargs)]
            append_runs_to_hdf5(filename, new_runs)
            h5open(filename, "r") do file
                @test keys(file) == ["configs", "initial_data", "run_1", "run_2"]
                @test length(keys(file["run_2"])) == 4
            end
        end

        # Test case 3: save_obj_to_hdf5
        @testset "save_obj_to_hdf5" begin
            filename = joinpath(dir, "test2.h5")
            obj = CustomType(6, 7.0, "save_obj")
            h5open(filename, "w") do file
                Skylight.save_obj_to_hdf5(file, "custom_type", obj)
                @test keys(file) == ["custom_type"]
                @test length(keys(file["custom_type"])) == 4
            end
        end

        # Test case 4: save_nested_dict_to_hdf5
        @testset "save_nested_dict_to_hdf5" begin
            filename = joinpath(dir, "test3.h5")
            nested_dict = Dict(:key1 => "value1", :key2 => Dict(:key3 => "value3"))
            h5open(filename, "w") do file
                Skylight.save_nested_dict_to_hdf5(file, nested_dict)
                @test keys(file) == ["key1", "key2"]
                @test keys(file["key2"]) == ["key3"]
            end
        end

        @testset "to_hdf5_compatible_dict and is_hdf5_supported_type tests" begin
            test_dict = Dict(:a => 1,
                :b => 2.0,
                :c => 3 + 4im,
                :d => "test",
                :e => sin,
                :f => MyType2(1, 2.0, 3 + 4im, "test", sin),
                :g => nothing)
            test_obj = MyType2(1, 2.0, 3 + 4im, "test", sin)

            @testset "is_hdf5_supported_type" begin
                @test Skylight.is_hdf5_supported_type(1)
                @test Skylight.is_hdf5_supported_type(1.0)
                @test Skylight.is_hdf5_supported_type(1 + 2im)
                @test Skylight.is_hdf5_supported_type("test")
                @test !Skylight.is_hdf5_supported_type(sin)
                @test !Skylight.is_hdf5_supported_type(test_obj)
                @test !Skylight.is_hdf5_supported_type(test_dict)
            end

            @testset "to_hdf5_compatible_dict for Dict" begin
                hdf5_dict = Skylight.to_hdf5_compatible_dict(test_dict)

                @test hdf5_dict[:a] == 1
                @test hdf5_dict[:b] == 2.0
                @test hdf5_dict[:c] == 3 + 4im
                @test hdf5_dict[:d] == "test"
                @test hdf5_dict[:e] == "sin"
                @test hdf5_dict[:f] isa Dict
                @test hdf5_dict[:g] == "Nothing"
            end

            @testset "to_hdf5_compatible_dict for custom type" begin
                hdf5_dict = Skylight.to_hdf5_compatible_dict(test_obj)

                @test hdf5_dict["a"] == 1
                @test hdf5_dict["b"] == 2.0
                @test hdf5_dict["c"] == 3 + 4im
                @test hdf5_dict["d"] == "test"
                @test hdf5_dict["e"] == "sin"
            end
        end

        @testset "to_hdf5_compatible_dict for DECallback" begin
            # Assuming 'callback' is an instance of VectorContinuousCallback
            hdf5_dict = Skylight.to_hdf5_compatible_dict(callback)

            # '_typename' should be the type of the callback
            @test hdf5_dict["_typename"] == string(typeof(callback))

            # Check other fields
            # `condition`, `affect!`, `affect_neg!`, `initialize`, `finalize` fields are functions
            @test hdf5_dict["condition"] == string(getfield(callback, :condition))
            @test hdf5_dict["affect!"] == string(getfield(callback, :affect!))
            @test hdf5_dict["affect_neg!"] == string(getfield(callback, :affect_neg!))
            @test hdf5_dict["initialize"] == string(getfield(callback, :initialize))
            @test hdf5_dict["finalize"] == string(getfield(callback, :finalize))

            # `len` field is a number
            @test hdf5_dict["len"] == getfield(callback, :len)

            # `idxs` field can be `nothing` or `AbstractArray`
            if getfield(callback, :idxs) === nothing
                @test hdf5_dict["idxs"] == "nothing"
            else
                @test hdf5_dict["idxs"] == getfield(callback, :idxs)
            end

            # `save_positions` field is a Tuple of Bools converted to BitVector
            expected_str = string(:(BitVector($([Bool(x) for x in getfield(callback,
                :save_positions)]))))
            @test hdf5_dict["save_positions"] == expected_str

            # `interp_points`, `abstol`, `reltol` fields are numbers
            @test hdf5_dict["interp_points"] == getfield(callback, :interp_points)
            @test hdf5_dict["abstol"] == getfield(callback, :abstol)
            @test hdf5_dict["reltol"] == getfield(callback, :reltol)

            # `repeat_nudge` field is a Rational
            expected_str = string(:(Rational($(numerator(getfield(callback,
                    :repeat_nudge))),
                $(denominator(getfield(callback, :repeat_nudge))))))
            @test hdf5_dict["repeat_nudge"] == expected_str
        end
    end
end
