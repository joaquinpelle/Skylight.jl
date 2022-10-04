module Skylight

using Parameters

include("coordinatesystems.jl")
include("spacetimes.jl")
include("emissionmodels.jl")
include("initialdataconfigurations.jl")

include("utils.jl")
#include("polarcap.jl")
include("initialdataOTE.jl")

end
