module Skylight

using Parameters

include("spacetimes.jl")
include("emissionmodels.jl")
include("imageplane.jl")
include("initialdataconfigurations.jl")

include("utils.jl")
#include("polarcap.jl")
include("initialdataOTE.jl")

end
