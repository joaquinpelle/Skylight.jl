"""
    static_slice(u::SVector{N, T}, ::Val{start}, ::Val{length}) where {N, T, start, length}

Extract a slice of length `length` from the input SVector, starting at index `start`.
Be careful since it doesn't have any bounds checking!
"""
@generated function static_slice(u::SVector{N, T}, ::Val{start}, ::Val{length}) where {N, T, start, length}
    # Create expressions to extract the desired elements from the input SVector
    elements = ntuple(i -> :(u[$(i + start - 1)]), length)  
    # Construct a new SVector from these elements
    return Expr(:call, SVector{length, T}, elements...)
end

"""
    compose_N_times(N, f, x)

Compose function `f` to the input `x` for `N` times.

# Arguments
- `N`: The number of times to apply the function `f`.
- `f`: The function to be applied.
- `x`: The input to the function `f`.

# Example
```julia
matrix = rand(Bool, 5, 5)
compose_N_times(3, detect_edges, matrix)
```
"""
compose_N_times(N, f, x) = _compose_N_times(Val(N), f, x)

@generated function _compose_N_times(::Val{N}, f, x) where {N}
    expr = :(f(x))
    for _ in 1:(N - 1)
        expr = :(f($expr))
    end
    return expr
end
