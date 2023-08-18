include("cache.jl")
include("views.jl")
include("finalpositions.jl")
include("normalization.jl")
include("redshift.jl")
include("intensities.jl")
include("fluxes.jl")
include("spectrum.jl")
include("lineemission.jl")

function postprocess_init(initial_data::AbstractMatrix, output_data::AbstractMatrix, ::VacuumOTEConfigurations, ::Nothing)
    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have eight components.")) 
    nrays = size(initial_data, 2)
    array = zeros(nrays)
    return array
end

function postprocess_init(initial_data::AbstractMatrix, output_data::AbstractMatrix, ::VacuumOTEConfigurations, observation_energies::AbstractVector)
    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    eight_components(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have eight components.")) 
    nrays = size(initial_data, 2)
    NE = length(observation_energies)
    array = zeros(NE, nrays)
    return array
end

function postprocess_init(initial_data::AbstractMatrix, output_data::AbstractMatrix, ::NonVacuumOTEConfigurations)
    same_size(initial_data, output_data) || throw(DimensionMismatch("The initial and output data must have the same size."))
    nrays = size(initial_data, 2)
    NE = (size(initial_data, 1) - 8)/2
    array = zeros(NE,nrays)
    return array
end

task_init(configurations::VacuumOTEConfigurations; kwargs...) = task_init(configurations.camera, configurations; kwargs...)

function task_init(::ImagePlane, configurations::VacuumOTEConfigurations)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    cache = postprocess_cache(configurations)
    return spacetime, model, coords_top, cache
end

function task_init(camera::PinholeCamera, configurations::VacuumOTEConfigurations; observer_four_velocity)
    spacetime = configurations.spacetime
    model = configurations.radiative_model
    coords_top = coordinates_topology(spacetime)
    cache = postprocess_cache(configurations)
    observer_metric!(cache, camera.position, spacetime)
    observer_four_velocity!(cache, observer_four_velocity) 
    return spacetime, model, coords_top, cache
end
