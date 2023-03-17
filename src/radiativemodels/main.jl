abstract type RadiativeModel end

abstract type SurfaceEmissionModel <: RadiativeModel end

abstract type OpaqueInteriorSurfaceEmissionModel <: SurfaceEmissionModel end

abstract type NeutronStarHotSpots <: OpaqueInteriorSurfaceEmissionModel end
abstract type BlackHoleAccretionDisk <: SurfaceEmissionModel end
abstract type BlackHoleCorona <: SurfaceEmissionModel end

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("blackholeaccretiondisk.jl")
include("bosonstaraccretiondisk.jl")
include("staracrosswormhole.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")
include("utils.jl")
