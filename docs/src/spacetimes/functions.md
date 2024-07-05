# Geometric functions

### Allocating methods
```@docs
metric(position::AbstractVector, spacetime::AbstractSpacetime)
metric_inverse(position::AbstractVector, spacetime::AbstractSpacetime)
volume_element(position::AbstractVector, spacetime::AbstractSpacetime)
christoffel(position::AbstractVector, spacetime::AbstractSpacetime)
```

### Non-allocating methods
For functions called within tight loops, we suggest to use the following non-allocating methods 
```@docs
metric!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractSpacetimeCache})
metric_inverse!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
volume_element(position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
christoffel!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractChristoffelCache})
```
