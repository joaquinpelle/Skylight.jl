"""
    StarAcrossWormhole <: AbstractSurfaceEmissionModel

Uniform temperature star accross a traversable wormhole ([`ChargedWormholeSpacetimeRegularCoordinates`](@ref)).

# Fields
- `l_center::Float64`: The center of the star in the regular radial coordinate. 
- `star_radius::Float64`: The radius of the star.

# Examples
```julia
star = StarAcrossWormhole(l_center=5.0, star_radius=1.0)
```
"""
@with_kw struct StarAcrossWormhole <: AbstractSurfaceEmissionModel
    l_center::Float64
    star_radius::Float64
    @assert star_radius>0.0 "star_radius must be positive"
end

opaque_interior_surface_trait(::StarAcrossWormhole) = IsOpaqueInteriorSurface()
stationarity(::StarAcrossWormhole) = IsStationary()