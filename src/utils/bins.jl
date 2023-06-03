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
        if bin_index == 1 || bin_index === nothing
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

function create_bins(; bin_size::Union{Number,Nothing}=nothing, num_bins::Union{Int,Nothing}=nothing, start::Number, stop::Number)
    if !(bin_size===nothing)
        return start:bin_size:stop
    elseif !(num_bins===nothing)
        return range(start, stop=stop, length=num_bins+1)
    else
        error("Must provide either bin_size or num_bins")
    end
end