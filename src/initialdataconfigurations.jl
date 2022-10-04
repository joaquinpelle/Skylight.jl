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
    initial_times::Vector{Float64}

end

my_zeros(configurations) = zeros(8, number_of_initial_conditions(configurations))

get_initial_times(configurations) = configurations.initial_times

function number_of_initial_conditions(configurations::OTEInitialDataConfigurations)
    
    number_of_nodes = number_of_nodes(configurations.image_plane) 
    number_of_times = length(configurations.initial_times)
    
    return number_of_nodes*number_of_times 
    
end

function number_of_initial_conditions(configurations::ETOInitialDataConfigurations)
    
    number_of_points = get_number_of_points(emission_model)
    number_of_times = length(configurations.initial_times)
    return number_of_points*number_of_packets_per_point*number_of_times 
    
end
