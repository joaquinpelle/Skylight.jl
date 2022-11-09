export OTEConfigurations
export ETOConfigurations

abstract type Configurations end

abstract type TransferScheme end

struct ObserverToEmitter <: TransferScheme end
struct EmitterToObserver <: TransferScheme end

@with_kw struct OTEConfigurations{S<:Spacetime, M<:EmissionModel} <: Configurations
    
    spacetime::S
    emission_model::M
    image_plane::ImagePlane
    initial_times::Vector{Float64} 

end

transfer_scheme(configurations::OTEConfigurations) = ObserverToEmitter()

@with_kw struct ETOConfigurations{S<:Spacetime, M<:EmissionModel} <: Configurations
    
    spacetime::S
    emission_model::M
    number_of_packets_per_point::Int64
    observer_distance::Float64

end

transfer_scheme(configurations::ETOConfigurations) = EmitterToObserver()

my_zeros(configurations) = zeros(8, number_of_initial_conditions(configurations))

get_initial_times(configurations::OTEConfigurations) = configurations.initial_times

get_initial_data_cache(configurations::OTEConfigurations) = OTEInitialDataCache()
get_initial_data_cache(configurations::ETOConfigurations) = ETOInitialDataCache()


function get_initial_positions(configurations::ETOConfigurations)
    
    times = zero_times(configurations)
    space_positions = get_space_positions(configurations)
    
    return eachcol([times'; space_positions])

end

function get_space_positions(configurations::ETOConfigurations)
    
    coord_system = coordinate_system_class(configurations.spacetime)
    space_positions = get_space_positions(configurations.emission_model, coord_system)

    return space_positions

end

function zero_times(configurations::ETOConfigurations)
    
    npoints = get_number_of_points(configurations.emission_model)
    return repeat([0.0],npoints)

end

function number_of_initial_conditions(configurations::OTEConfigurations)
     
    number_of_times = length(configurations.initial_times)
    
    return number_of_nodes(configurations.image_plane)*number_of_times 
    
end

function number_of_initial_conditions(configurations::ETOConfigurations)
    
    number_of_points = get_number_of_points(configurations.emission_model)
    number_of_packets_per_point = configurations.number_of_packets_per_point

    return number_of_points*number_of_packets_per_point
    
end
