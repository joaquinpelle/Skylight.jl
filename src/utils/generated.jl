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
