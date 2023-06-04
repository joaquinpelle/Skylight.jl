include("emittertoobserver/cache.jl")
include("observertoemitter/cache.jl")
include("emittertoobserver/emittertoobserver.jl")
include("observertoemitter/imageplane.jl")
include("observertoemitter/pinholecamera.jl")

get_initial_data(configurations::VacuumOTEConfigurations) = get_initial_data(configurations.camera, configurations)