module Skylight

using Parameters

include("geometry.jl")
include("types.jl")

include("spacetimes/minkowski.jl")
include("spacetimes/schwarzschild.jl")
include("spacetimes/kerr.jl")
include("spacetimes/johannsen.jl")
include("spacetimes/charged_wormhole.jl")

include("initialdataOTE.jl")

end
