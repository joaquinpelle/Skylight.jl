 """
    save_to_hdf5(filename, configurations, initial_data, runs)

Save the provided data to an HDF5 file. The data includes configurations, initial_data, and a collection of runs. Each run contains output_data, callback, callback_parameters, and kwargs for the solver function.

If the file exists, it uses the "r+" mode to open the file for reading and writing. If the file does not exist, it uses the "w" mode to create a new file.

The file will be organized as follows:
- /configurations: A group containing the configurations data.
- /initial_data: A dataset containing the initial data.
- /runs: A group containing subgroups for each run. Each subgroup will have:
    - /output_data: A dataset containing the output data.
    - /callback: A subgroup containing the callback custom type data.
    - /callback_parameters: A subgroup containing the callback parameters data.
    - /kwargs: A subgroup containing the kwargs data.

# Arguments
* `filename`: The name of the HDF5 file to save the data to. If the file already exists, new data will be appended without overwriting existing content.
* `configurations`: A custom type containing the configurations data. It will be converted to a dictionary using the `to_hdf5_compatible_dict` function.
* `initial_data`: The initial data for the equations.
* `runs`: An array of tuples containing output_data, callback custom type, callback_parameters, and kwargs for each run. The callback custom type and callback_parameters should be custom types, which will be converted to dictionaries using the `to_hdf5_compatible_dict` function.

# Returns
* Nothing.
"""

function save_to_hdf5(filename::String, configurations::AbstractConfigurations, initial_data::Array, runs::Array{T,}) where {T<:Run}
    
    mode = isfile(filename) ? "r+" : "w"

    h5open(filename, mode) do file

        # Save initial data
        write(file, "initial_data", initial_data)
        # Save configurations
        save_obj_to_hdf5(file, "configs", configurations)
        # Save runs
        save_runs_to_hdf5(file, runs)
    end
end

"""
    append_runs_to_hdf5(filename::String, runs::Vector{Tuple{AbstractArray, Any, Any, Dict{Symbol,}}})

Appends a set of runs to an existing HDF5 file containing simulation data. The runs are saved
under separate groups named "run_N", where N is the run index. The input `runs` is a vector of
tuples, each containing the output data, callback, callback parameters, and keyword arguments
(kwargs) for a single run.

# Arguments
- `filename::String`: The path to the HDF5 file where the runs will be appended.
- `runs::Vector{Tuple{AbstractArray, Any, Any, Dict{Symbol,}}}`: A vector of tuples containing
  the output data, callback, callback parameters, and kwargs for each run to be appended.

# Example
```julia
filename = "simulation_data.h5"
runs = [
    (output_data1, callback1, callback_parameters1, kwargs1),
    (output_data2, callback2, callback_parameters2, kwargs2),
]
append_runs_to_hdf5(filename, runs)
"""

function append_runs_to_hdf5(filename::String, runs::Array{T,}) where {T<:Run}
    
    h5open(filename, "r+") do file
        # Count the number of existing runs
        num_runs = count(k -> startswith(k, "run_"), keys(file))

        # Save runs
        save_runs_to_hdf5(file, runs, num_runs=num_runs)

    end
end

"""
    save_runs_to_hdf5(file::HDF5File, runs::Vector{Tuple{AbstractArray, Any, Any, Dict{Symbol,}}}; num_runs::Int=0)

Saves a set of runs to an open HDF5 file containing simulation data. The runs are saved
under separate groups named "run_N", where N is the run index. The input `runs` is a vector of
tuples, each containing the output data, callback, callback parameters, and keyword arguments
(kwargs) for a single run. The optional `num_runs` argument specifies the starting index for
the run numbering.

# Arguments
- `file::HDF5File`: An open HDF5 file where the runs will be saved.
- `runs::Vector{Tuple{AbstractArray, Any, Any, Dict{Symbol,}}}`: A vector of tuples containing
  the output data, callback, callback parameters, and kwargs for each run to be saved.
- `num_runs::Int` (optional, default: 0): The starting index for the run numbering.

# Example
```julia
using HDF5

filename = "simulation_data.h5"
runs = [
    (output_data1, callback1, callback_parameters1, kwargs1),
    (output_data2, callback2, callback_parameters2, kwargs2),
]

h5open(filename, "w") do file
    save_runs_to_hdf5(file, runs)
end
"""

function save_runs_to_hdf5(file, runs; num_runs=0)
    
    for (i, run) in enumerate(runs)

        run_group = create_group(file, "run_$(num_runs + i)")
        
        # Unpack the run tuple
        output_data, callback, callback_parameters, kwargs = to_tuple(run)
        
        # Save output data, callback, callback_parameters, and kwargs to the new run group
        write(run_group, "output_data", output_data)
        save_obj_to_hdf5(run_group, "callback", callback)
        save_obj_to_hdf5(run_group, "callback_parameters", callback_parameters)
        save_obj_to_hdf5(run_group, "kwargs", kwargs)
    end

end

"""
    save_obj_to_hdf5(group::HDF5File, name::String, obj::Any)

Save an object to an HDF5 group by converting it to an HDF5 compatible dictionary and then saving the dictionary to the group. The function creates a new subgroup with the given name and saves the converted dictionary using the `save_nested_dict_to_hdf5` function.

# Arguments
* `group`: An HDF5 group to which the object should be saved.
* `name`: A string representing the name of the subgroup to create for the object.
* `obj`: An object of any type that should be saved to the HDF5 group.
"""
function save_obj_to_hdf5(group, name, obj)
    
    subgroup = create_group(group, name)
    dict = to_hdf5_compatible_dict(obj)
    save_nested_dict_to_hdf5(subgroup, dict)

end

"""
    save_nested_dict_to_hdf5(group::HDF5File, nested_dict::Dict{Symbol,})

Recursively save a nested dictionary to an HDF5 group. The function iterates through the key-value pairs in the input dictionary. If a value is a dictionary, it creates a new subgroup and recursively saves the nested dictionary. If a value is not a dictionary, it is saved directly to the current group.

# Arguments
* `group`: An HDF5 group to which the nested dictionary should be saved.
* `nested_dict`: A nested dictionary with keys of type `Symbol` and values of any type.
"""
function save_nested_dict_to_hdf5(group, nested_dict)
    for (key, value) in nested_dict
        if isa(value, Dict)
            subgroup = create_group(group, string(key))
            save_nested_dict_to_hdf5(subgroup, value)
        else
            group[string(key)] = value
        end
    end
end

"""
    to_hdf5_compatible_dict(dict::Dict{Symbol,}, max_depth::Int) -> Dict{Symbol,}

Recursively convert a Julia dictionary into an HDF5 compatible dictionary, up to a specified maximum depth. The function checks each value in the input dictionary to see if it is an HDF5 supported type, and if not, it attempts to convert the value into a compatible format.

The resulting HDF5 compatible dictionary can be used to save data to an HDF5 file using the HDF5.jl package.

# Arguments
* `dict`: A Julia dictionary with keys of type `Symbol` and values of any type.
* `max_depth`: The maximum depth allowed for recursion when converting nested dictionaries.

# Returns
* An HDF5 compatible dictionary with keys of type `Symbol` and values of supported types.
"""

function to_hdf5_compatible_dict(dict::Dict{T, S}; depth::Int=0, max_depth::Int=8 ) where {T, S}
    if depth > max_depth
        return nothing
    end

    hdf5_dict = Dict{T, S}()

    for (key, value) in dict
        if isa(value, NoSaveField)
            continue
        elseif is_hdf5_supported_type(value)
            hdf5_dict[key] = value
        elseif isa(value, Dict)
            hdf5_dict[key] = to_hdf5_compatible_dict(value, depth=depth+1, max_depth=max_depth)
        elseif isa(value, Function)
            hdf5_dict[key] = string(value)
        elseif value === nothing
            hdf5_dict[key] = "Nothing"
        else
            fields_dict = to_hdf5_compatible_dict(value, depth=depth+1, max_depth=max_depth)
            if !isempty(fields_dict)
                hdf5_dict[key] = fields_dict
            end
        end
    end

    return hdf5_dict
end

"""
    to_hdf5_compatible_dict(obj::T; depth::Int=0, max_depth::Int=5) where T

Recursively convert a custom type `obj` to a dictionary. This function is useful for preparing custom types for saving to an HDF5 file when the types are not directly supported by the HDF5.jl package.

The function will attempt to convert each field of the custom type to a dictionary if the field is not an HDF5-supported type. The conversion process is limited by a maximum recursion depth, specified by the `max_depth` parameter, which defaults to 5.

If a field is nothing, the value will be saved as the string "nothing".

If a field is a function, the function will be converted to a string using the `string` function.

If a field is not an HDF5-supported type and has no fields that can be converted, the function will return `nothing` when the maximum recursion depth is reached.

# Arguments
* `obj::T`: The custom type instance to be converted to a dictionary.
* `depth::Int=0` (optional): The current recursion depth, starting at 0. This parameter is used internally and should not typically be set by the user.
* `max_depth::Int=5` (optional): The maximum recursion depth allowed for converting nested custom types.

# Returns
* A `Dict{String, Any}` containing the fields of the input custom type, with nested custom types recursively converted to dictionaries. The dictionary will also include a "_typename" key with the name of the custom type as a string.
"""
function to_hdf5_compatible_dict(obj::T; depth::Int=0, max_depth::Int=8) where T
    if depth > max_depth
        return nothing
    end

    d = Dict{String, Any}()
    d["_typename"] = string(T)
    for field in fieldnames(T)
        value = getfield(obj, field)
        
        # Check if the field is no-save
        if isa(value, NoSaveField)
            continue
        # Check if the field is of type Nothing
        elseif isa(value, Nothing)
            value = "nothing"
        # Check if the field is a subtype of Function
        elseif isa(value, Function)
            value = string(value)
        # Check if the field is a custom type not supported by HDF5.jl
        elseif !is_hdf5_supported_type(value)
            value = to_hdf5_compatible_dict(value, depth=depth+1, max_depth=max_depth)
        end
        
        d[string(field)] = value
    end
    return d
end

"""
    to_hdf5_compatible_dict(cb::T; depth::Int=0, max_depth::Int=8) where {T<:SciMLBase.DECallback}

Converts a DifferentialEquations.jl callback into a dictionary that can be saved to an HDF5 file
for later instantiation.

# Arguments
- `cb`: The callback to be converted into a dictionary.
- `depth`: The current depth of recursion. This is used to prevent infinite recursion. 
  (default: 0)
- `max_depth`: The maximum allowed depth of recursion. If the `depth` exceeds this value, 
  the function returns `nothing`. (default: 8)

# Returns
- A dictionary where the keys are the names of the fields of the object and the values are 
  the corresponding field values. Field values are recursively converted to dictionaries 
  if they are of custom types not supported by HDF5.jl.

# Notes
- The dictionary includes an additional entry with the key "_typename" and the name of the type 
  of the object as the value. This is used to reconstruct the original object when loading 
  the data from the HDF5 file.
- The function specifically handles certain fields of the DECallback type.
"""
function to_hdf5_compatible_dict(cb::T; depth::Int=0, max_depth::Int=8) where {T<:SciMLBase.DECallback}
    if depth > max_depth
        return nothing
    end

    d = Dict{String, Any}()
    d["_typename"] = string(T)
    for field in fieldnames(T)
        value = getfield(cb, field)
        
        if field == :repeat_nudge
            value = string(:(Rational($(numerator(value)), $(denominator(value)))))
        elseif field == :save_positions
            value = string(:(BitVector($([Bool(x) for x in value])))) 
        # Check if the field is of type Nothing
        elseif isa(value, Nothing)
            value = "nothing"
        # Check if the field is a subtype of Function
        elseif isa(value, Function)
            value = string(value)
        # Check if the field is a custom type not supported by HDF5.jl
        elseif !is_hdf5_supported_type(value)
            value = to_hdf5_compatible_dict(value, depth=depth+1, max_depth=max_depth)
        end
        
        d[string(field)] = value
    end
    return d
end

"""
    is_hdf5_supported_type(value)

Check if a given value is of a type that can be directly saved using the HDF5.jl package. This function can be used to determine whether a custom type or one of its fields can be saved to an HDF5 file without further conversion.

According to the HDF5.jl documentation, the supported types include:
- Signed and unsigned integers (8, 16, 32, and 64 bits)
- Float32 and Float64
- Complex versions of numeric types
- Arrays of supported numeric types (including complex versions)
- AbstractString (ASCIIString and UTF8String)
- Arrays of supported string types (ASCIIString and UTF8String)

# Arguments
* `value`: The value to be checked for HDF5 support.

# Returns
* A `Bool` indicating whether the input value is of a type supported by HDF5.jl.
"""
function is_hdf5_supported_type(value)
    supported_types = [Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32, Float64]
    supported_complex_types = [Complex{t} for t in supported_types]
    supported_string_types = [String]
    all_supported_types = vcat(supported_types, supported_complex_types, supported_string_types)

    T = isa(value, Array) ? eltype(value) : typeof(value)
    
    return T in all_supported_types
end