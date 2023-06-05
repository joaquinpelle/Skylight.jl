@with_kw struct StarAcrossWormhole <: AbstractSurfaceEmissionModel 

    l_center::Float64
    star_radius::Float64
    @assert star_radius > 0.0 "star_radius must be positive"
end

opaque_interior_surface_trait(::StarAcrossWormhole) = IsOpaqueInteriorSurface()
