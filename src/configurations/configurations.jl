export NonVacuumOTEConfigurations
export VacuumOTEConfigurations
export VacuumETOConfigurations

abstract type Configurations end

abstract type OTEConfigurations <: Configurations end
abstract type ETOConfigurations <: Configurations end

@with_kw struct NonVacuumOTEConfigurations{S<:Spacetime, M<:RadiativeModel} <: OTEConfigurations
    
    spacetime::S
    radiative_model::M
    image_plane::ImagePlane
    observed_times::Vector{Float64}
    observed_energies::Vector{Float64}
    unit_mass_in_solar_masses::Float64

end

@with_kw struct VacuumOTEConfigurations{S<:Spacetime, M<:RadiativeModel} <: OTEConfigurations
    
    spacetime::S
    radiative_model::M
    image_plane::ImagePlane
    observed_times::Vector{Float64}
    unit_mass_in_solar_masses::Float64

end

@with_kw struct VacuumETOConfigurations{S<:Spacetime, M<:RadiativeModel} <: ETOConfigurations
    
    spacetime::S
    radiative_model::M
    number_of_points::Int64
    number_of_packets_per_point::Int64
    observer_distance::Float64
    unit_mass_in_solar_masses::Float64

end

VacuumConfigurations = Union{VacuumOTEConfigurations, VacuumETOConfigurations}
NonVacuumConfigurations = Union{NonVacuumOTEConfigurations,}

function my_zeros(configurations::NonVacuumConfigurations)

    NE = length(configurations.observed_energies)
    
    return zeros(8+2*NE, number_of_initial_conditions(configurations))

end

my_zeros(configurations::VacuumConfigurations) = zeros(8, number_of_initial_conditions(configurations))

get_observed_times(configurations::OTEConfigurations) = configurations.observed_times

get_initial_data_cache(configurations::OTEConfigurations) = OTEInitialDataCache()
get_initial_data_cache(configurations::ETOConfigurations) = ETOInitialDataCache()

get_postprocess_cache(configurations::OTEConfigurations) = OTEPostProcessCache()

function get_initial_positions(configurations::ETOConfigurations)
    
    times = zero_times(configurations)
    space_positions = get_space_positions(configurations)
    
    return eachcol([times'; space_positions])

end

function get_space_positions(configurations::ETOConfigurations)
    
    npoints = configurations.number_of_points

    coord_system = coordinate_system_class(configurations.spacetime)
    space_positions = get_space_positions(npoints, configurations.radiative_model, coord_system)

    return space_positions

end

function zero_times(configurations::ETOConfigurations)
    
    npoints = configurations.number_of_points
    return repeat([0.0],npoints)

end

function number_of_initial_conditions(configurations::OTEConfigurations)
     
    number_of_times = length(configurations.observed_times)
    
    return number_of_nodes(configurations.image_plane)*number_of_times 
    
end

function number_of_initial_conditions(configurations::ETOConfigurations)
    
    number_of_points = configurations.number_of_points
    number_of_packets_per_point = configurations.number_of_packets_per_point

    return number_of_points*number_of_packets_per_point
    
end
