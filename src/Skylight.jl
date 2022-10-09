module Skylight

using LinearAlgebra
using Parameters

include("spacetimes/main.jl")
include("emissionmodels/main.jl")
include("imageplane.jl")
include("initialdataconfigurations.jl")

include("utils/main.jl")
include("initialdataOTE.jl")

end
