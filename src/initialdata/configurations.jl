export OTEInitialDataConfigurations
export ETOInitialDataConfigurations

abstract type InitialDataConfigurations end

@with_kw struct OTEInitialDataConfigurations{S<:Spacetime} <: InitialDataConfigurations
    
    spacetime::S
    image_plane::ImagePlane 
    initial_times::Vector{Float64}

end

@with_kw struct ETOInitialDataConfigurations{S<:Spacetime, M<:EmissionModel} <: InitialDataConfigurations
    
    spacetime::S
    emission_model::M
    number_of_packets_per_point::Int64

end

my_zeros(configurations) = zeros(8, number_of_initial_conditions(configurations))

get_initial_times(configurations::OTEInitialDataConfigurations) = configurations.initial_times

get_cache(configurations::OTEInitialDataConfigurations) = OTEInitialDataCache()
get_cache(configurations::ETOInitialDataConfigurations) = ETOInitialDataCache()


function get_initial_positions(configurations::ETOInitialDataConfigurations)
    
    times = zero_times(configurations)
    space_positions = get_space_positions(configurations)
    
    return eachcol([times'; space_positions])

end

function get_space_positions(configurations::ETOInitialDataConfigurations)
    
    coord_system = coordinate_system_kind(configurations.spacetime)
    space_positions = get_space_positions(configurations.emission_model, coord_system)

    return space_positions

end

function zero_times(configurations::ETOInitialDataConfigurations)
    
    npoints = get_number_of_points(configurations.emission_model)
    return repeat([0.0],npoints)

end


function number_of_initial_conditions(configurations::OTEInitialDataConfigurations)
     
    number_of_times = length(configurations.initial_times)
    
    return number_of_nodes(configurations.image_plane)*number_of_times 
    
end

function number_of_initial_conditions(configurations::ETOInitialDataConfigurations)
    
    number_of_points = get_number_of_points(configurations.emission_model)
    number_of_packets_per_point = configurations.number_of_packets_per_point

    return number_of_points*number_of_packets_per_point
    
end
