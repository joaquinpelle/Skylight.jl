export ImagePlane
export OTEInitialDataConfigurations

abstract type InitialDataConfigurations end

@with_kw struct ImagePlane

    observer_distance :: Float64
    observer_inclination_in_degrees :: Float64
    horizontal_side_image_plane :: Float64
    vertical_side_image_plane :: Float64
    horizontal_number_of_nodes :: Int32
    vertical_number_of_nodes :: Int32
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)

end

@with_kw struct OTEInitialDataConfigurations{S<:Spacetime} <: InitialDataConfigurations
    
    spacetime::S
    image_plane::ImagePlane 
    initial_times::Vector{Float64}

end

@with_kw struct ETOInitialDataConfigurations{S<:Spacetime, M<:EmissionModel} <: InitialDataConfigurations
    
    spacetime::S
    emission_model::M
    packets_per_point::Int64
    initial_times::Vector{Float64}

end

get_initial_times(configurations) = configurations.initial_times
