include("macros/types.jl")
include("spacetimes/types.jl")
include("radiativemodels/types.jl")
include("configurations/types.jl")
include("initialdata/types.jl")
include("transfer/types.jl")
include("postprocess/types.jl")


SkylightCustomType = Union{AbstractSpacetime, AbstractRadiativeModel, AbastractCallbackParameters}

