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
    create_bins(; bin_size::Number=NaN, num_bins::Int=NaN, start_range::Number=NaN, 
                 end_range::Number=NaN)

Create bins for binning data.

# Keywords
- `bin_size::Number=NaN`: Size of each bin. Either `bin_size` or `num_bins` must be specified.
- `num_bins::Int=NaN`: Number of bins. Either `bin_size` or `num_bins` must be specified.
- `start_range::Number`: Lower bound of the range to be binned. 
- `end_range::Number`: Upper bound of the range to be binned. 

# Returns
- `bins`: Array of the bin edges.
"""
    
function create_bins(; bin_size::Number=NaN, num_bins::Int=NaN, start_range::Number, end_range::Number)
    if !isnan(bin_size)
        return start_range:bin_size:end_range
    elseif !isnan(num_bins)
        return range(start_range, stop=end_range, length=num_bins+1)
    else
        error("Must provide either bin_size or num_bins")
    end
end