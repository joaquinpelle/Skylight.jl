abstract type EmissionModel end

abstract type SurfaceEmissionModel <: EmissionModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end
abstract type TransparentSurfaceEmissionModel <: SurfaceEmissionModel end

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("thinaccretiondisk.jl")


function set_unit_surface_normal!(vector, position, gcache, spacetime, model, coord_system)

    set_surface_differential!(vector, position, model, coord_system)

    set_metric_inverse!(gcache, position, spacetime)
    raise_index!(vector,gcache)

    set_metric!(gcache,position,spacetime)
    normalize_spacelike!(vector, gcache)

end
