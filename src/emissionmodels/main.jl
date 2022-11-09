abstract type EmissionModel end

abstract type SurfaceEmissionModel <: EmissionModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end

abstract type NeutronStarHotSpots <: OpaqueInteriorSurfaceEmissionModel end
abstract type BlackHoleAccretionDisk <: SurfaceEmissionModel end
abstract type BlackHoleCorona <: SurfaceEmissionModel end

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("blackholeaccretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("utils.jl")
