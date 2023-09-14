# Geometric functions

```@docs
metric(position::AbstractVector, spacetime::AbstractSpacetime)
metric_inverse(position::AbstractVector, spacetime::AbstractSpacetime)
volume_element(position::AbstractVector, spacetime::AbstractSpacetime)
christoffel(position::AbstractVector, spacetime::AbstractSpacetime)
```

### Non-allocating methods
For non-allocating versions to be used within tight loops, use the following methods 
```@docs
metric!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractSpacetimeCache})
metric_inverse!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
volume_element(position::AbstractVector, spacetime::AbstractSpacetime, g::AbstractMatrix, cache::Union{Nothing, AbstractSpacetimeCache})
christoffel!(metric::AbstractMatrix, position::AbstractVector, spacetime::AbstractSpacetime, cache::Union{Nothing, AbstractChristoffelCache})
```

#### Cache allocation

```@docs
allocate_cache(spacetime::AbstractSpacetime)
allocate_christoffel_cache(spacetime::AbstractSpacetime)
```