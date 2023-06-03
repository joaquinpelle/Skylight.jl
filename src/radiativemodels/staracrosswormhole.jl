@with_kw struct StarAcrossWormhole <: AbstractSurfaceEmissionModel 

    l_center::Float64
    star_radius::Float64

end

opaque_interior_surface_trait(::StarAcrossWormhole) = IsOpaqueInteriorSurface()
