abstract type EmissionModel end

abstract type SurfaceEmissionModel <: EmissionModel end
abstract type VolumeEmissionModel <: EmissionModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end
abstract type TransparentSurfaceEmissionModel <: SurfaceEmissionModel end

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("thinaccretiondisk.jl")


function set_unit_surface_normal!(vector, position, metric, model, coord_system)

    set_surface_differential!(vector, position, model, coord_system)
    vector .= raise_index(vector,metric)
    normalize_spacelike!(vector, metric)

end
