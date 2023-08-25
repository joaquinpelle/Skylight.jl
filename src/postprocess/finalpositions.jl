"""
    is_final_position_at_source(output_data, configurations)

Determine whether the final position of each ray is at the source for a given set of configurations.

# Arguments
- `output_data`: a multi-dimensional array containing output data for each ray.
- `configurations`: a data structure containing the spacetime and radiative model configurations.

# Returns
A boolean grid indicating whether each ray's final position is at the source.

# Notes
The function `is_final_position_at_source` should already be defined elsewhere, and is used within this function to check the final position of each ray.
"""
function is_final_position_at_source(output_data::AbstractMatrix, configurations)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    nrays = size(output_data, 2)
    at_source = zeros(Bool, nrays)
    @threads for i in axes(output_data, 2)
        @views pf = output_data[1:4, i]
        at_source[i] = is_final_position_at_source(pf, spacetime, model)
    end
    return grid_view(at_source, configurations)
end

"""
    is_final_position_at_edge(output_data, configurations)

Determine whether the final position of each ray is at the edge of the source for a given set of configurations.

# Arguments
- `output_data`: a multi-dimensional array containing output data for each ray.
- `configurations`: a data structure containing the spacetime and radiative model configurations.

# Returns
A boolean grid indicating whether each ray's final position is at the edge of the source.
"""
function is_final_position_at_edge(output_data, configurations)
    at_source = is_final_position_at_source(output_data, configurations)
    return detect_edges(at_source)
end
"""
    is_final_position_at_edge(width::Int, output_data, configurations)

Determine whether the final position of each ray is at the edge of the source for a given set of configurations, considering a specified edge width.

# Arguments
- `width`: the width (in number of cells) to consider as the edge of the source.
- `output_data`: a multi-dimensional array containing output data for each ray.
- `configurations`: a data structure containing the spacetime and radiative model configurations.

# Returns
A boolean grid indicating whether each ray's final position is within the specified width of the edge of the source.
"""
function is_final_position_at_edge(width::Int, output_data, configurations)
    at_source = is_final_position_at_source(output_data, configurations)
    return detect_edges(width, at_source)
end

"""
    is_final_position_at_observer(output_data, configurations)

Determine whether the final position of each ray is at the observer for a given set of configurations.

# Arguments
- `output_data`: a multi-dimensional array containing output data for each ray.
- `configurations`: a data structure containing the spacetime and radiative model configurations.

# Returns
A boolean grid indicating whether each ray's final position is at the observer.

# Notes
The function `is_final_position_at_observer` should already be defined elsewhere, and is used within this function to check the final position of each ray.
"""
function is_final_position_at_observer(output_data, configurations)
    nrays = size(output_data, 2)
    at_observer = zeros(Bool, nrays)
    @threads for i in axes(output_data, 2)
        @views pf = output_data[1:4, i]
        at_observer[i] = is_final_position_at_observer(pf, configurations)
    end
    return at_observer
end
