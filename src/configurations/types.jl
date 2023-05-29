abstract type AbstractConfigurations end
abstract type AbstractOTEConfigurations <: AbstractConfigurations end
abstract type AbstractETOConfigurations <: AbstractConfigurations end

abstract type AbstractCamera end

@with_kw struct ImagePlane <: AbstractCamera
    distance::Float64
    observer_inclination_in_degrees::Float64
    horizontal_side::Float64
    vertical_side::Float64
    horizontal_number_of_nodes::Int
    vertical_number_of_nodes::Int
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)
end

@with_kw struct PinholeCamera <: AbstractCamera 
    position::Vector{Float64}
    four_velocity::Vector{Float64}
    horizontal_aperture::Float64 #This is distance*cos(horizontal_aperture_angle)
    vertical_aperture::Float64   #This is distance*cos(horizontal_aperture_angle)
    horizontal_number_of_nodes::Int
    vertical_number_of_nodes::Int
    direction::Vector{Float64} = to_center(position, four_velocity)
end

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