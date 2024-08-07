include("camera.jl")

function isvacuum(::AbstractConfigurations)
    error("isvacuum not implemented for this type of configurations")
end
isvacuum(::VacuumOTEConfigurations) = Vacuum()
isvacuum(::VacuumETOConfigurations) = Vacuum()
isvacuum(::NonVacuumOTEConfigurations) = NonVacuum()

stationarity(configurations::AbstractConfigurations) = stationarity(configurations.spacetime) && stationarity(configurations.radiative_model)
spherical_symmetry(configurations::AbstractConfigurations) = spherical_symmetry(configurations.spacetime) && spherical_symmetry(configurations.radiative_model)
axisymmetry(configurations::AbstractConfigurations) = axisymmetry(configurations.spacetime) && axisymmetry(configurations.radiative_model)
helical_symmetry(configurations::AbstractConfigurations) = helical_symmetry(configurations.spacetime) && helical_symmetry(configurations.radiative_model)

function my_zeros(configurations::AbstractConfigurations)
    my_zeros(isvacuum(configurations), configurations)
end
my_zeros(::Vacuum, configurations) = zeros(8, number_of_initial_conditions(configurations))
function my_zeros(::NonVacuum, configurations)
    NE = length(configurations.observation_energies)
    return zeros(8 + 2 * NE, number_of_initial_conditions(configurations))
end

axes_ranges(configurations::AbstractOTEConfigurations) = axes_ranges(configurations.camera)

function initial_data_cache(configurations::AbstractOTEConfigurations)
    initial_data_cache(configurations.camera, configurations.spacetime)
end
function initial_data_cache(configurations::AbstractETOConfigurations)
    ETOInitialDataCache(configurations.spacetime, configurations.radiative_model)
end
function postprocess_cache(configurations::AbstractOTEConfigurations)
    postprocess_cache(configurations.camera,
        configurations.spacetime,
        configurations.radiative_model)
end

function postprocess_cache(configurations::AbstractETOConfigurations)
    spacetime_cache = allocate_cache(configurations.spacetime)
    model_cache = allocate_cache(configurations.radiative_model)
    return ETOPostProcessCache(spacetime_cache=spacetime_cache, model_cache=model_cache)
end

function initial_positions(configurations::AbstractETOConfigurations)
    times = zero_times(configurations)
    return eachcol([times'; space_positions(configurations)])
end

function space_positions(configurations::AbstractETOConfigurations)
    npoints = configurations.number_of_points
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    scache = allocate_cache(spacetime)
    coords_top = coordinates_topology(configurations.spacetime)
    return space_positions(npoints, spacetime, model, coords_top, scache)
end

function zero_times(configurations::AbstractETOConfigurations)
    npoints = configurations.number_of_points
    return repeat([0.0], npoints)
end

function number_of_initial_conditions(configurations::AbstractOTEConfigurations)
    number_of_initial_conditions(configurations.camera)
end

function number_of_initial_conditions(configurations::AbstractETOConfigurations)
    number_of_points = configurations.number_of_points
    number_of_packets_per_point = configurations.number_of_packets_per_point
    return number_of_points * number_of_packets_per_point
end

function max_radius(configurations::AbstractOTEConfigurations)
    max_radius(configurations.camera, configurations.spacetime)
end
max_radius(configurations::AbstractETOConfigurations) = configurations.max_radius
