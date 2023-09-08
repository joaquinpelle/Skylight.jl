"""
    bin_values_and_sum_weights(bins, values, weights)

Bin `values` and sum `weights` in each bin.

# Arguments
- `bins`: Array of bin edges.
- `values`: Array of values to be binned.
- `weights`: Array of weights to be summed in each bin.

# Returns
- `binned_values`: Array of the sum of weights in each bin.

# Notes
- The length of `values` and `weights` must be the same.
- Values outside the range of `bins` are ignored.
"""
function bin_values_and_sum_weights(bins, values, weights)
    if length(values) != length(weights)
        throw(ArgumentError("Length of values and weights must be the same"))
    end

    binned_values = zeros(length(bins) - 1)

    for i in eachindex(values)
        value = values[i]
        weight = weights[i]

        # find the bin index for the current value
        bin_index = findfirst(b -> b > value, bins)

        # skip values outside the bin range
        if bin_index == 1 || isnothing(bin_index)
            continue
        end

        # decrement by 1 to get the correct index
        bin_index -= 1

        # increment the sum for this bin
        binned_values[bin_index] += weight
    end

    return binned_values
end

"""
    average_inside_bins(q, x, bins)

Calculate the averages of a quantity `q` at value `x` in each bin defined by `bins`.
"""
function average_inside_bins(q, x, bins)
    same_size(q, x) || throw(ArgumentError("q and x must have the same size"))
    # Initialize an array to hold the sum of `q` values in each bin
    qsums = zeros(eltype(q), length(bins)-1)
    # Initialize an array to hold the count of `q` values in each bin
    qcounts = zeros(Int, length(bins)-1)
    averages = zero(qsums)
    for (xvalue, qvalue) in zip(x, q)
        # Find the bin index for the current radii value
        bin_index = searchsortedlast(bins, xvalue)

        # If the bin_index is valid (i.e., the radius value is not beyond the last bin edge)
        if bin_index > 0 && bin_index < length(bins)
            # Update the sum and count of `q` values for this bin
            qsums[bin_index] += qvalue
            qcounts[bin_index] += 1
        end
    end
    # Calculate the average `q` values for each bin
    for i in 1:(length(bins)-1)
        if qcounts[i] == 0
            continue
        end
        averages[i] = qsums[i] / qcounts[i]
    end
    return averages
end

"""
    create_bins(; bin_size::Union{Number, Nothing}=nothing, num_bins::Union{Int, Nothing}=nothing, 
                 start::Number, stop::Number)

Create bins for binning data.

# Keywords
- `bin_size::Union{Number, Nothing}=nothing`: Size of each bin. Either `bin_size` or `num_bins` must be specified.
- `num_bins::Union{Int, Nothing}=nothing`: Number of bins. Either `bin_size` or `num_bins` must be specified.
- `start::Number`: Lower bound of the range to be binned. 
- `stop::Number`: Upper bound of the range to be binned. 

# Returns
- `bins`: Array of the bin edges.

# Notes
- If `bin_size` is provided, the function will create bins with that size from `start` to `stop`.
- If `num_bins` is provided, the function will create that number of equally spaced bins from `start` to `stop`.
- If neither `bin_size` nor `num_bins` is provided, the function will throw an error.
"""

function create_bins(;
    bin_size::Union{Number, Nothing} = nothing,
    num_bins::Union{Int, Nothing} = nothing,
    start::Number,
    stop::Number)
    if !(bin_size === nothing)
        return start:bin_size:stop
    elseif !(num_bins === nothing)
        return range(start, stop = stop, length = num_bins + 1)
    else
        error("Must provide either bin_size or num_bins")
    end
end

"""
    infer_num_bins(q, at_source, start, stop, bin_size_conditioner, edge_width, camera)

Infer the number of bins from the energy quotients and the bin size conditioner. It chooses the larges number of bins that satisfies
that the bin size is larger than the conditioner times the maximum local variation of the energy quotients, for a given approximation of the local variation.

# Arguments

- `q`: Array of energy quotients.
- `at_source`: Array of booleans indicating whether the final position is at the source.
- `start`: Lower bound of the range to be binned.
- `stop`: Upper bound of the range to be binned.
- `bin_size_conditioner`: Conditioner of the bin size.
- `edge_width`: Width of the edge to be ignored for bin size conditioning (since edges have unusually large local variations, especially in the presence of light-rings and similar phenomena).
- `camera`: Image plane.

# Returns

- `num_bins`: Number of bins.
"""
function infer_num_bins(q, at_source, start, stop, bin_size_conditioner, edge_width, camera)
    Nα, Nβ = numbers_of_pixels_per_side(camera)
    dα, dβ = grid_spacing(camera)

    at_edge = detect_edges(edge_width, reshape(at_source, Nα, Nβ))

    δq = approximate_gradient_norm(reshape(q, Nα, Nβ), dα, dβ)
    max_δq = maximum(δq[.!at_edge])
    return Int(floor((stop - start) / (bin_size_conditioner * max_δq)))
end

"""
    midpoints(edges)

Returns midpoints of given edges
"""
midpoints(edges::AbstractVector) = 0.5*(edges[1:end-1] + edges[2:end])
widths(edges::AbstractVector) = edges[2:end] - edges[1:end-1]