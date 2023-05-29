abstract type AbstractConfigurations end
abstract type AbstractOTEConfigurations <: AbstractConfigurations end
abstract type AbstractETOConfigurations <: AbstractConfigurations end


@with_kw struct NonVacuumOTEConfigurations{S<:AbstractSpacetime, M<:AbstractRadiativeModel} <: AbstractOTEConfigurations
    
    spacetime::S
    radiative_model::M
    image_plane::ImagePlane
    observed_times::Vector{Float64}
    observed_energies::Vector{Float64}
    unit_mass_in_solar_masses::Float64

end

@with_kw struct VacuumOTEConfigurations{S<:AbstractSpacetime, M<:AbstractRadiativeModel} <: AbstractOTEConfigurations
    
    spacetime::S
    radiative_model::M
    image_plane::ImagePlane
    observed_times::Vector{Float64}
    unit_mass_in_solar_masses::Float64

end

@with_kw struct VacuumETOConfigurations{S<:AbstractSpacetime, M<:AbstractRadiativeModel} <: AbstractETOConfigurations
    
    spacetime::S
    radiative_model::M
    number_of_points::Int
    number_of_packets_per_point::Int
    observer_distance::Float64
    unit_mass_in_solar_masses::Float64

end

VacuumConfigurations = Union{VacuumETOConfigurations, VacuumOTEConfigurations}
NonVacuumConfigurations = Union{NonVacuumOTEConfigurations,}