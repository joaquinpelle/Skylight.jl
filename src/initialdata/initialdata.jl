include("emittertoobserver/cache.jl")
include("observertoemitter/cache.jl")
include("emittertoobserver/emittertoobserver.jl")
include("observertoemitter/imageplane.jl")
include("observertoemitter/pinholecamera.jl")

get_initial_data(configurations::VacuumOTEConfigurations) = get_initial_data(configurations.camera, configurations)

function dump_metric_in!(cache::AbstractInitialDataCache, position, spacetime::AbstractSpacetime)
    set_metric!(cache.metric, position, spacetime)
end