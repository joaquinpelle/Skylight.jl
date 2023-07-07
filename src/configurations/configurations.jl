include("camera.jl")

isvacuum(::AbstractConfigurations) = error("isvacuum not implemented for this type of configurations")
isvacuum(::VacuumOTEConfigurations) = Vacuum()
isvacuum(::VacuumETOConfigurations) = Vacuum()
isvacuum(::NonVacuumOTEConfigurations) = NonVacuum()

my_zeros(configurations::AbstractConfigurations) = my_zeros(isvacuum(configurations), configurations)
my_zeros(::Vacuum, configurations) = zeros(8, number_of_initial_conditions(configurations))
function my_zeros(::NonVacuum, configurations)
    NE = length(configurations.observation_energies)
    return zeros(8+2*NE, number_of_initial_conditions(configurations))
end

axes_ranges(configurations::AbstractOTEConfigurations) = axes_ranges(configurations.camera)

initial_data_cache(::AbstractOTEConfigurations) = initial_data_cache(configurations.camera, configurations.spacetime)
initial_data_cache(configurations::AbstractETOConfigurations) = ETOInitialDataCache(configurations.spacetime, configurations.radiative_model)
postprocess_cache(configurations::AbstractOTEConfigurations) = postprocess_cache(configurations.camera, configurations.spacetime, configurations.model)

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

number_of_initial_conditions(configurations::AbstractOTEConfigurations) = number_of_initial_conditions(configurations.camera)

function number_of_initial_conditions(configurations::AbstractETOConfigurations)
    number_of_points = configurations.number_of_points
    number_of_packets_per_point = configurations.number_of_packets_per_point
    return number_of_points*number_of_packets_per_point
end

max_radius(configurations::AbstractOTEConfigurations) = max_radius(configurations.camera, configurations.spacetime) 
max_radius(configurations::AbstractETOConfigurations) = configurations.observer_distance