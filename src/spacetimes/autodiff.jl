""""
Calculates the Christoffel symbols of a given spacetime metric using the forward-mode automatic differentiation package ForwardDiff.

Parameters:
- Γ₂: mutable array of size (4,4,4) to store the resulting Christoffel symbols.
- position: tuple of four numbers representing a position in spacetime.
- spacetime: object representing the spacetime.
- cache: object of type AutoDiffChristoffelCache with containers for metric elements and derivatives.

Returns: nothing.

Note: for the automatic differentiation to work on a given spacetime, any cache array stored in the spacetime
struct for metric calculations must be wrapped by the DiffCache() method as:
    
    ```
    @with_kw struct KerrSpacetimeKerrSchildCoordinates{T} <: AbstractSpacetime

        M::Float64
        a::Float64
    
        @assert M >= 0.0
        @assert abs(a) <= M 
    
        #Metric cache
        l::T = DiffCache(zeros(4))
    
    end
    ```

This is because automatic differentiation keeps two versions of each variable, a Real and a Dual version, the latter being used
to compute derivatives at each node of the chain rule. Also, in the `metric!(g, position, spacetime)` function, the caches must be accessed via the function `get_tmp` as in
`get_tmp(spacetime.l, position)`, so that the appropriate version of the cache is returned according to the element type of position
when `metric!` is called.

"""
function christoffel!(Γ₂, position, spacetime::AbstractSpacetime, cache::AutoDiffChristoffelCache)
    
    g = cache.g
    ginv = cache.ginv
    ∂g = cache.∂g
    spacetime_metric_field = cache.spacetime_metric_field
    cfg = cache.cfg

    metric_inverse!(ginv, position, spacetime, g)
    metric_jacobian!(∂g, position, spacetime_metric_field, g, cfg)
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
    metric_jacobian!(∂g, position, spacetime::AbstractSpacetime, g)

Computes the Jacobian matrix of the metric function with respect to spacetime coordinates using forward-mode automatic differentiation.

Parameters:
- ∂g: mutable array of size (4,4,4) to store the resulting Jacobian matrix.
- position: tuple of four numbers representing a position in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) containing the metric evaluated at the given position.

Returns: nothing.
"""
function metric_jacobian!(∂g, position, spacetime::AbstractSpacetime, g)
    ForwardDiff.jacobian!(∂g, metric_field(spacetime), g, position)
    return nothing
end

"""
    metric_jacobian!(∂g, position, spacetime_metric_field::Function, g, cfg::ForwardDiff.JacobianConfig)

Computes the Jacobian matrix of the metric function with respect to spacetime coordinates using forward-mode automatic differentiation.

Parameters:
- ∂g: mutable array of size (4,4,4) to store the resulting Jacobian matrix.
- position: tuple of four numbers representing a position in spacetime.
- spacetime: object representing the spacetime.
- g: array of size (4,4) containing the metric evaluated at the given position.
- cfg: object of type ForwardDiff.JacobianConfig with preallocated work buffers 

Returns: nothing.
"""
function metric_jacobian!(∂g, position, spacetime_metric_field::Function, g, cfg::ForwardDiff.JacobianConfig)
    ForwardDiff.jacobian!(∂g, spacetime_metric_field, g, position, cfg)
    return nothing
end