include("emittertoobserver/cache.jl")
include("observertoemitter/cache.jl")
include("emittertoobserver/emittertoobserver.jl")
include("observertoemitter/imageplane.jl")
include("observertoemitter/pinholecamera.jl")

function initialize(configurations::AbstractOTEConfigurations)
    initialize(configurations.camera, configurations)
end
