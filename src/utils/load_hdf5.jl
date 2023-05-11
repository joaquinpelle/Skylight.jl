"""
    load_initial_data_from_hdf5(filename::String)

Load initial data from an HDF5 file.

# Arguments
- `filename::String`: The path to the HDF5 file.

# Returns
- The initial data stored in the HDF5 file.
"""
function load_initial_data_from_hdf5(filename::String)
    h5open(filename, "r") do file
        return read(file, "initial_data")
    end
end

"""
    load_configurations_from_hdf5(filename::String)

Load configurations data from an HDF5 file.

# Arguments
- `filename::String`: The path to the HDF5 file.

# Returns
- An instance of the configuration type containing the configurations data stored in the HDF5 file.
"""
function load_configurations_from_hdf5(filename::String)
    h5open(filename, "r") do file
        configs_group = file["configs"]
        configs_dict = load_nested_dict_from_hdf5(configs_group)

        return instantiate_custom_type(configs_dict)
    end
end

"""
    load_runs(filename::String, run_indices::Vector{Int})

Load a set of runs specified by their indices from an HDF5 file.

# Arguments
- `filename`: The path to the HDF5 file.
- `run_indices`: A vector containing the indices of the runs to be loaded.

# Returns
- A vector of tuples, each containing the output_data, callback dictionary,
  callback_parameters, and kwargs dictionary for the specified runs.
"""
function load_runs_from_hdf5(filename::String, run_indices::Vector{Int})
    h5open(filename, "r") do file
        runs = []
        for i in run_indices
            run_group = file["run_$(i)"]
            output_data = read(run_group, "output_data")
            
            cb = load_callback_from_hdf5(filename, i)
            cbp = load_callback_params_from_hdf5(filename, i)
            
            kwargs_group = file["run_$(run_index)/kwargs"]
            kwargs_dict = load_nested_dict_from_hdf5(kwargs_group)
            
            push!(runs, (output_data, cb, cbp, kwargs_dict))
        end
        return runs
    end
end

"""
    load_output_data_from_hdf5(filename::String, run_indices::Vector{Int})

Load output data for a selected set of runs from an HDF5 file.

# Arguments
- `filename::String`: The path to the HDF5 file.
- `run_indices::Vector{Int}`: A vector of integers indicating which runs' output data to load.

# Returns
- An array containing the output data for the selected runs.
"""
function load_output_data_from_hdf5(filename::String, run_indices::Vector{Int})
    h5open(filename, "r") do file
        runs_data = []
        for i in run_indices
            run_group = file["run_$(i)"]
            output_data = read(run_group, "output_data")
            push!(runs_data, output_data)
        end
        return runs_data
    end
end

"""
    load_output_data_from_hdf5(filename::String, run_index::Int)

Load the output data of a specific run from an HDF5 file.

# Arguments
- `filename`: The path to the HDF5 file.
- `run_index`: The index of the run to load the output data from.

# Returns
- The output data of the specified run.
"""
function load_output_data_from_hdf5(filename::String, run_index::Int)
    h5open(filename, "r") do file
        run_group = file["run_$(run_index)"]
        output_data = read(run_group, "output_data")
        return output_data
    end
end

"""
    load_callback_params_from_hdf5(filename::String, run_index::Int)

Load the callback parameters of a specific run from an HDF5 file and
instantiate the custom type if possible.

# Arguments
- `filename`: The path to the HDF5 file.
- `run_index`: The index of the run to load the callback parameters from.

# Returns
- An instance of the custom type representing the callback parameters, or
  a dictionary if the custom type cannot be instantiated.
"""
function load_callback_params_from_hdf5(filename::String, run_index::Int)
    h5open(filename, "r") do file
        cbp_group = file["run_$(run_index)/callback_parameters"]
        cbp_dict = load_nested_dict_from_hdf5(cbp_group)

        return instantiate_custom_type(cbp_dict)
    end
end
"""
    load_callback_from_hdf5(filename::String, run_index::Int)

Load the callback dictionary of a specific run from an HDF5 file.

# Arguments
- `filename`: The path to the HDF5 file.
- `run_index`: The index of the run to load the callback from.

# Returns
- A dictionary representing the callback of the specified run.
"""

function load_callback_from_hdf5(filename::String, run_index::Int)
    h5open(filename, "r") do file
        cb_group = file["run_$(run_index)/callback"]
        return load_nested_dict_from_hdf5(cb_group)
    end
end


"""
    load_kwargs_from_hdf5(filename::String, run_index::Int)

Load the kwargs dictionary of a specific run from an HDF5 file.

# Arguments
- `filename`: The path to the HDF5 file.
- `run_index`: The index of the run to load the kwargs from.

# Returns
- A dictionary representing the kwargs of the specified run.
"""

function load_kwargs_from_hdf5(filename::String, run_index::Int)
    h5open(filename, "r") do file
        kwargs_group = file["run_$(run_index)/kwargs"]
        return load_nested_dict_from_hdf5(kwargs_group)
    end
end
"""
    instantiate_custom_type(dict)

Create an instance of a custom type using a dictionary. This version supports custom types
created using the `@with_kw` macro from the Parameters.jl package and handles nested
dictionaries. The dictionary and its nested subdictionaries must contain a `_typename` key whose value is the name of the custom type. 

# Arguments
- `dict::Dict{Symbol, Any}`: A dictionary containing the custom type's name (stored in the `_typename` key) and field values.

# Returns
- An instance of the custom type with the field values specified in the input dictionary, including instances created from nested dictionaries.

# Example

```julia
using Parameters

@with_kw struct MyNestedType
    a::Int
end

@with_kw struct MyTypeWithKW
    b::Float64
    nested::MyNestedType
end

my_dict = Dict(:_typename => "MyTypeWithKW",
               :b => 2.0,
               :nested => Dict(:_typename => "MyNestedType",
                               :a => 1))

my_instance = instantiate_custom_type(my_dict)
"""
function instantiate_custom_type(dict::Dict{Symbol, Any})
    typename = dict[:_typename]
    T = eval(Meta.parse(typename))

    # Instantiate nested custom types for any fields that are also dictionaries
    for (k, v) in dict
        if isa(v, Dict)
            dict[k] = instantiate_custom_type(v)
        end
    end

    kwarg_dict = Dict(Symbol(k) => v for (k, v) in pairs(dict) if (k != :_typename))
    
    return T(; kwarg_dict...)
end

function instantiate_custom_type(dict::Dict{String, Any})
    # Convert the dictionary keys to symbols
    symbol_dict = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in dict)
    return instantiate_custom_type(symbol_dict)
end

"""
    load_nested_dict_from_hdf5(group::HDF5Group)

Convert an HDF5 group to a dictionary.

# Arguments
- `group::HDF5Group`: The HDF5 group to convert to a dictionary.

# Returns
- A dictionary containing the data stored in the HDF5 group.
"""
function load_nested_dict_from_hdf5(group::HDF5.Group)
    nested_dict = Dict{Symbol, Any}()
    
    for name in keys(group)
        obj = group[name]

        if isa(obj, HDF5.Group)
            nested_dict[Symbol(name)] = load_nested_dict_from_hdf5(obj)
        elseif isa(obj, HDF5.Dataset)
            nested_dict[Symbol(name)] = read(obj)
        else
            @warn "Unsupported object type found in HDF5 group: $name"
        end
    end
    
    return nested_dict
end
