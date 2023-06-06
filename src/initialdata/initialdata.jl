include("emittertoobserver/cache.jl")
include("observertoemitter/cache.jl")
include("emittertoobserver/emittertoobserver.jl")
include("observertoemitter/imageplane.jl")
include("observertoemitter/pinholecamera.jl")

initialize(configurations::VacuumOTEConfigurations) = initialize(configurations.camera, configurations)