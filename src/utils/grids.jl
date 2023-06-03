"""
    approximate_gradient_norm(z::Array, dx_spacing::Real, dy_spacing::Real)

Compute the norm of the gradient of a function given a 2D grid of function values `z` and grid spacings `dx_spacing` and `dy_spacing`.

The gradient is computed using central differences in the interior and forward or backward differences at the edges.

# Arguments
- `z::Array`: 2D grid of function values.
- `dx_spacing::Real`: Grid spacing in the x direction.
- `dy_spacing::Real`: Grid spacing in the y direction.

# Returns
- `grad_norm::Array`: 2D grid of gradient magnitudes.

"""
function approximate_gradient_norm(z::Array, dx_spacing::Real, dy_spacing::Real)
    dx = (circshift(z, (-1, 0)) - circshift(z, (1, 0))) / (2*dx_spacing)
    dy = (circshift(z, (0, -1)) - circshift(z, (0, 1))) / (2*dy_spacing)

    dx[1, :] = (z[2, :] - z[1, :]) / dx_spacing
    dx[end, :] = (z[end, :] - z[end-1, :]) / dx_spacing
    dy[:, 1] = (z[:, 2] - z[:, 1]) / dy_spacing
    dy[:, end] = (z[:, end] - z[:, end-1]) / dy_spacing

    grad_norm = sqrt.(dx.^2 + dy.^2)
    return grad_norm
end

"""
    detect_edges(arr::Array{Bool, 2})

Detects edges in a two-dimensional Boolean array by comparing each element with its
neighbors in all four cardinal directions (up, down, left, and right). An edge is
considered to be present at a location if the Boolean value at that location is different
from the value of any of its neighbors.

# Arguments
- `arr::Array{Bool, 2}`: a two-dimensional Boolean array.

# Returns
- `Array{Bool, 2}`: a two-dimensional Boolean array of the same size as the input, where `true` indicates the presence of an edge and `false` indicates its absence.

# Notes
- This function assumes that the input array is padded with `false` (or 0) at the boundaries. You might want to add boundary condition checking or padding if your input does not meet this condition.
- This function uses "edge wrapping" when it shifts the array, meaning the values at the edge of the array "wrap around" to the opposite edge. If you want a different boundary condition (such as padding with `false`), you would need to adjust the creation of the shifted arrays accordingly.

"""
function detect_edges(arr::Union{Matrix{Bool}, BitMatrix})
    rows, cols = size(arr)

    # Create shifted arrays
    arr_shifted_up = vcat(arr[2:end, :], fill(false, 1, cols))
    arr_shifted_down = vcat(fill(false, 1, cols), arr[1:end-1, :])
    arr_shifted_left = hcat(arr[:, 2:end], fill(false, rows, 1))
    arr_shifted_right = hcat(fill(false, rows, 1), arr[:, 1:end-1])

    # Compute the sum of the squares of differences
    diff = (arr .- arr_shifted_up).^2 + 
           (arr .- arr_shifted_down).^2 + 
           (arr .- arr_shifted_left).^2 + 
           (arr .- arr_shifted_right).^2

    # Return the edges (non-zero elements in the diff array)
    return diff .> 0
end

"""
    detect_edges(N, arr)

Detects values that are at distance at most N from an edge. Composes `detect_edges` to 
the input `arr` for 1 to `N` times, and combines the results using an elementwise OR operation. 
Each application of `detect_edges` is performed on the result of the previous application and the
results are combined in succession.

# Arguments
- `N`::Int: The maximum number of times to apply the function `f`.
- `arr`::Union{Matrix{Bool}, BitMatrix}: The boolean array.

# Example
```julia
matrix = rand(Bool, 5, 5)
detect_edges(3, matrix)
````
In this example, the detect_edges function is applied 1 to 3 times to a 5x5 random
Boolean matrix, and the results are combined with an OR operation.
"""
detect_edges(N::Int, arr::Union{Matrix{Bool}, BitMatrix}) = _detect_edges(Val(N), arr)

@generated function _detect_edges(::Val{N}, arr) where N
    result = :(detect_edges(arr))
    for _ in 2:N
        result = :($result .| detect_edges($result))
    end
    return result
end

