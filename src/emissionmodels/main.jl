abstract type EmissionModel end

abstract type SurfaceEmissionModel <: EmissionModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("thinaccretiondisk.jl")
include("blackholecorona.jl")
include("starbehindwormhole.jl")
include("dummyextendedregion.jl")


function set_unit_surface_normal!(vector, position, metric, metric_inverse, model, coord_system)

    set_surface_differential!(vector, position, model, coord_system)
    vector .= raise_index(vector,metric_inverse)
    normalize_spacelike!(vector, metric)

end
