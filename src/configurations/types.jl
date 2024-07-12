abstract type AbstractConfigurations end
abstract type AbstractOTEConfigurations <: AbstractConfigurations end
abstract type AbstractETOConfigurations <: AbstractConfigurations end

abstract type AbstractCamera end

"""
    ImagePlane <: AbstractCamera

A type representing an image plane camera for ray-tracing simulations. This type assumes that light rays arrive parallel from a distant source.

# Fields

- `distance::Float64`: The distance of the observation point from the origin.
- `observer_inclination_in_degrees::Float64`: The inclination of the observer with respect to the z-axis in degrees.
- `horizontal_side::Float64`: The horizontal side length of the image plane.
- `vertical_side::Float64`: The vertical side length of the image plane.
- `horizontal_number_of_pixels::Int`: The number of pixels along the horizontal side of the image plane.
- `vertical_number_of_pixels::Int`: The number of pixels along the vertical side of the image plane.
- `observation_times::Vector{Float64}`: The observation times, defaulting to `[0.0]`.

# Constructor

```julia
camera = ImagePlane(distance = 500.0,
    observer_inclination_in_degrees = 90,
    horizontal_side = 23.0,
    vertical_side = 23.0,
    horizontal_number_of_pixels = 100,
    vertical_number_of_pixels = 100)
```
"""
@with_kw struct ImagePlane <: AbstractCamera
    distance::Float64
    observer_inclination_in_degrees::Float64
    horizontal_side::Float64
    vertical_side::Float64
    horizontal_number_of_pixels::Int
    vertical_number_of_pixels::Int
    observation_times::Vector{Float64} = [0.0]
    @assert 0 < distance
    @assert 0 <= observer_inclination_in_degrees <= 180
    @assert 0 < horizontal_side
    @assert 0 < vertical_side
    @assert 0 < horizontal_number_of_pixels
    @assert 0 < vertical_number_of_pixels
    observer_inclination_in_radians::Float64 = deg2rad(observer_inclination_in_degrees)
end

"""
    PinholeCamera <: AbstractCamera

A type representing a pinhole camera for ray-tracing simulations. This type allows for arbitrary positions and velocities of the observation point, providing a more general approach compared to the image plane.

# Fields

- `position::Vector{Float64}`: The position of the camera in spacetime.
- `horizontal_aperture_in_degrees::Float64`: The horizontal aperture of the camera in degrees.
- `vertical_aperture_in_degrees::Float64`: The vertical aperture of the camera in degrees.
- `horizontal_number_of_pixels::Int`: The number of pixels along the horizontal aperture.
- `vertical_number_of_pixels::Int`: The number of pixels along the vertical aperture.
- `four_velocity::Vector{Float64}`: The four-velocity of the observation frame, defaulting to the static frame once the initial data is created.

# Constructor

```julia
camera = PinholeCamera(position = [0.0, distance, π / 2 - π / 20, 0.0],
    horizontal_aperture_in_degrees = rad2deg(315 / distance),
    vertical_aperture_in_degrees = rad2deg(315 / distance),
    horizontal_number_of_pixels = 600,
    vertical_number_of_pixels = 600)
```
"""
@with_kw struct PinholeCamera <: AbstractCamera
    position::Vector{Float64}
    horizontal_aperture_in_degrees::Float64
    vertical_aperture_in_degrees::Float64
    horizontal_number_of_pixels::Int
    vertical_number_of_pixels::Int
    four_velocity::Vector{Float64} = zeros(4)
    @assert 0 < horizontal_aperture_in_degrees <= 180
    @assert 0 < vertical_aperture_in_degrees <= 90
    @assert 0 < horizontal_number_of_pixels
    @assert 0 < vertical_number_of_pixels
    horizontal_aperture_in_radians::Float64 = deg2rad(horizontal_aperture_in_degrees)
    vertical_aperture_in_radians::Float64 = deg2rad(vertical_aperture_in_degrees)
end

"""
    NonVacuumOTEConfigurations{S,M,C} <: AbstractOTEConfigurations

Cconfigurations object for non-vacuum observer-to-emitter (OTE) radiative transport problems.

# Fields
- `spacetime::S`: The spacetime in which the transport problem is set. Must be a subtype of `AbstractSpacetime`.
- `radiative_model::M`: The radiative model describing the source of radiation. Must be a subtype of `AbstractRadiativeModel` and non-vacuum.
- `camera::C`: The observational setup, including position and orientation. Must be a subtype of `AbstractCamera`.
- `observation_energies::Vector{Float64}`: A vector of energies (in ergs) at which the specific intensity is to be computed. All energies must be positive.
- `unit_mass_in_solar_masses::Float64`: The unit mass in solar masses used to determine the geometrized code units. Must be positive.

# Constructor
```julia
configurations = NonVacuumOTEConfigurations(
    spacetime = spacetime,
    radiative_model = model,
    camera = camera,
    observation_energies = [1e-8, 1e-7],
    unit_mass_in_solar_masses = 1.0)
```
"""
@with_kw struct NonVacuumOTEConfigurations{
    S <: AbstractSpacetime,
    M <: AbstractRadiativeModel,
    C <: AbstractCamera,
} <: AbstractOTEConfigurations
    spacetime::S
    radiative_model::M
    camera::C
    observation_energies::Vector{Float64}
    unit_mass_in_solar_masses::Float64
    @assert all(0 .< observation_energies) "all observation energies must be positive"
    @assert 0<unit_mass_in_solar_masses "unit_mass_in_solar_masses must be positive"
    @assert isa(isvacuum(radiative_model), NonVacuum) "radiative_model must be non-vacuum"
end

"""
    VacuumOTEConfigurations{S,M,C} <: AbstractOTEConfigurations

Cconfigurations object for vacuum observer-to-emitter (OTE) radiative transport problems.

# Fields
- `spacetime::S`: The spacetime in which the transport problem is set. Must be a subtype of `AbstractSpacetime`.
- `radiative_model::M`: The radiative model describing the source of radiation. Must be a subtype of `AbstractRadiativeModel` and non-vacuum.
- `camera::C`: The observational setup, including position and orientation. Must be a subtype of `AbstractCamera`.
- `unit_mass_in_solar_masses::Float64`: The unit mass in solar masses used to determine the geometrized code units. Must be positive.

# Constructor
```julia
configurations = NonVacuumOTEConfigurations(
    spacetime = spacetime,
    radiative_model = model,
    camera = camera,
    unit_mass_in_solar_masses = 1.0)
```
"""
@with_kw struct VacuumOTEConfigurations{
    S <: AbstractSpacetime,
    M <: AbstractRadiativeModel,
    C <: AbstractCamera,
} <: AbstractOTEConfigurations
    spacetime::S
    radiative_model::M
    camera::C
    unit_mass_in_solar_masses::Float64
    @assert 0<unit_mass_in_solar_masses "unit_mass_in_solar_masses must be positive"
    @assert isa(isvacuum(radiative_model), Vacuum) "radiative_model must be vacuum"
end

@with_kw struct VacuumETOConfigurations{
    S <: AbstractSpacetime,
    M <: AbstractRadiativeModel,
} <: AbstractETOConfigurations
    spacetime::S
    radiative_model::M
    number_of_points::Int
    number_of_packets_per_point::Int
    max_radius::Float64
    unit_mass_in_solar_masses::Float64
    @assert 0<number_of_points "number_of_points must be positive"
    @assert 0<number_of_packets_per_point "number_of_packets_per_point must be positive"
    @assert 0<max_radius "max_radius must be positive"
    @assert 0<unit_mass_in_solar_masses "unit_mass_in_solar_masses must be positive"
    @assert isa(isvacuum(radiative_model), Vacuum) "radiative_model must be vacuum"
end

abstract type AbstractTransferMethod end
struct ObserverToEmitter <: AbstractTransferMethod end
struct EmitterToObserver <: AbstractTransferMethod end
