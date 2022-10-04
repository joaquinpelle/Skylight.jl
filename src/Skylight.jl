module Skylight

using Parameters

include("abstracttypes.jl")
include("coordinatesystems.jl")
include("spacetimes.jl")
include("utils.jl")

include("emissionmodels.jl")
include("polarcap.jl")

include("initialdataOTE.jl")

end
