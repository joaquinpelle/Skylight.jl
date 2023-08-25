struct BogdanovPolarCap <: AbstractSurfaceEmissionModel end

opaque_interior_surface_trait(::BogdanovPolarCap) = IsOpaqueInteriorSurface()
