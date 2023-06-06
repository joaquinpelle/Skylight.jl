include("camera.jl")

function my_zeros(configurations::NonVacuumConfigurations)
    NE = length(configurations.observation_energies)
    return zeros(8+2*NE, number_of_initial_conditions(configurations))
end

my_zeros(configurations::VacuumConfigurations) = zeros(8, number_of_initial_conditions(configurations))

observation_times(configurations::AbstractOTEConfigurations) = configurations.camera.observation_times

pixel_coordinates_vectors(configurations::AbstractOTEConfigurations) = pixel_coordinates_vectors(configurations.camera)

initial_data_cache(configurations::AbstractOTEConfigurations) = initial_data_cache(configurations.camera)
initial_data_cache(::AbstractETOConfigurations) = ETOInitialDataCache()
postprocess_cache(configurations::AbstractOTEConfigurations) = postprocess_cache(configurations.camera)

function initial_positions(configurations::AbstractETOConfigurations)
    times = zero_times(configurations)
    return eachcol([times';  space_positions(configurations)])
end

function space_positions(configurations::AbstractETOConfigurations)
    npoints = configurations.number_of_points
    coords_top = coordinates_topology(configurations.spacetime)
    return space_positions(npoints, configurations.radiative_model, coords_top)
end

function zero_times(configurations::AbstractETOConfigurations)
    npoints = configurations.number_of_points
    return repeat([0.0],npoints)
end

function number_of_initial_conditions(configurations::AbstractOTEConfigurations)
    number_of_times = length(configurations.observation_times)
    return number_of_pixels(configurations.camera)*number_of_times  
end

function number_of_initial_conditions(configurations::AbstractETOConfigurations)
    number_of_points = configurations.number_of_points
    number_of_packets_per_point = configurations.number_of_packets_per_point
    return number_of_points*number_of_packets_per_point
end

max_radius(configurations::AbstractOTEConfigurations) = max_radius(configurations.camera, configurations.spacetime) 
max_radius(configurations::AbstractETOConfigurations) = configurations.observer_distance