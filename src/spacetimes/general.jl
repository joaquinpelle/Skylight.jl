"""
Calculates the Christoffel symbols of a given spacetime metric using the forward-mode automatic differentiation package ForwardDiff.

Parameters:
- Γ₂: mutable array of size (4,4,4) to store the resulting Christoffel symbols.
- position: tuple of four numbers representing a point in spacetime.
- spacetime: object representing the spacetime.
- cache: object of type GeneralChristoffelCache with containers for metric elements and derivatives.

Returns: nothing.

Note: for the automatic differentiation to work on a given spacetime, any cache array stored in the spacetime
struct for metric calculations must be wrapped by the DiffCache() method as:
    
    ```
    @with_kw struct KerrSpacetimeKerrSchildCoordinates{T} <: BlackHoleSpacetime

        M::Float64
        a::Float64
    
        @assert M >= 0.0
        @assert abs(a) <= M 
    
        #Metric cache
        l::T = DiffCache(zeros(4))
    
    end
    ```

This is because automatic differentiation keeps two versions of each variable, a Real and a Dual version, the latter being used
to compute derivatives at each node of the chain rule. Also, in the `set_metric!(g, position, spacetime)` function, the caches must be accessed via the function `get_tmp` as in
`get_tmp(spacetime.l, position)`, so that the appropriate version of the cache is returned according to the element type of position
when `set_metric!` is called.

"""
function set_christoffel!(Γ₂, position, spacetime, cache::GeneralChristoffelCache)
    
    g = cache.g
    ginv = cache.ginv
    ∂g = cache.∂g

    set_metric_inverse!(ginv, position, spacetime, g)
    set_metric_jacobian!(∂g, position, spacetime, g)
    @inbounds begin
        for k in 1:4
            for j in 1:4
                for i in 1:j
                    for l in 1:4
                        Γ₂[l, i, j] += 0.5 * ginv[k, l] * (∂g[k, j, i] + ∂g[i, k, j] - ∂g[i, j, k])
                    end
                end
            end
        end
        for j in 1:4
            for i in j+1:4
                for l in 1:4
                    Γ₂[l, i, j] = Γ₂[l, j, i]
                end
            end
        end
    end
    return nothing
end

"""
Computes the Jacobian matrix of the metric function with respect to spacetime coordinates using forward-mode automatic differentiation.

Parameters:
- ∂g: mutable array of size (4,4,4) to store the resulting Jacobian matrix.
- position: tuple of four numbers representing a point in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) containing the metric evaluated at the given position.

Returns: nothing.
"""
function set_metric_jacobian!(∂g, position, spacetime, g)
    reshape(ForwardDiff.jacobian!(∂g, (g,q) -> set_metric!(g,q,spacetime), g, position), 4, 4, 4)
    return nothing
end

"""
Computes the inverse of the given metric at the given position using a fast inversion
for 4x4 symmetric matrices.

Parameters:
- ginv: mutable array of size (4,4) to store the resulting inverse metric.
- position: tuple of four numbers representing a point in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) to store the metric evaluated at the given position.

Returns: nothing.
"""
function set_metric_inverse!(ginv, position, spacetime, g)
    set_metric!(g, position, spacetime)
    inverse_4x4_symmetric!(ginv, g)
    return nothing
end

@with_kw mutable struct GeneralChristoffelCache
    g::Array{Float64,2} = zeros(4,4)
    ginv::Array{Float64,2} = zeros(4,4)
    ∂g::Array{Float64,3} = zeros(4,4,4)
end

allocate_christoffel_cache(spacetime) = GeneralChristoffelCache()