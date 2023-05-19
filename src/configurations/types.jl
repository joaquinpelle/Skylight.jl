@with_kw struct ImagePlane

    distance :: Float64
    observer_inclination_in_degrees :: Float64
    horizontal_side_image_plane :: Float64
    vertical_side_image_plane :: Float64
    horizontal_number_of_nodes :: Int64
    vertical_number_of_nodes :: Int64
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)

end

abstract type AbstractConfigurations end

abstract type OTEConfigurations <: AbstractConfigurations end
abstract type ETOConfigurations <: AbstractConfigurations end

@with_kw struct NonVacuumOTEConfigurations{S, M} <: OTEConfigurations
    
    spacetime::S
    radiative_model::M
    image_plane::ImagePlane
    observed_times::Vector{Float64}
    observed_energies::Vector{Float64}
    unit_mass_in_solar_masses::Float64

end

@with_kw struct VacuumOTEConfigurations{S, M} <: OTEConfigurations
    
    spacetime::S
    radiative_model::M
    image_plane::ImagePlane
    observed_times::Vector{Float64}
    unit_mass_in_solar_masses::Float64

end

@with_kw struct VacuumETOConfigurations{S, M} <: ETOConfigurations
    
    spacetime::S
    radiative_model::M
    number_of_points::Int64
    number_of_packets_per_point::Int64
    observer_distance::Float64
    unit_mass_in_solar_masses::Float64

end

VacuumConfigurations = Union{VacuumOTEConfigurations, VacuumETOConfigurations}
NonVacuumConfigurations = Union{NonVacuumOTEConfigurations,}