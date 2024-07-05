#Required
function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    model,
    coords_top,
    spacetime_cache,
    model_cache)
    error("rest_frame_four_velocity! not defined for this model.")
end
function rest_frame_four_velocity!(vector, position, metric, spacetime, model, coords_top)
    error("rest_frame_four_velocity! not defined for this model.")
end
function rest_frame_bolometric_intensity(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top,
    cache)
    error("rest_frame_bolometric_intensity for this model.")
end
function rest_frame_bolometric_intensity(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top)
    error("rest_frame_bolometric_intensity for this model.")
end
function rest_frame_specific_intensity(position,
    momentum,
    energy,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top,
    cache)
    error("rest_frame_specific_intensity not defined for this model.")
end
function rest_frame_specific_intensity(position,
    momentum,
    energy,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top)
    error("rest_frame_specific_intensity not defined for this model.")
end

"""
    line_emission_profile(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)

    Emissivity profile for line emission radiative models.

#Arguments

- `position::`: Position of the emission.
- `momentum::`: Momentum of the emission (frequency and direction). Must be null.
- `rest_frame_four_velocity::`: Rest frame four velocity of the model at `position`. Must but be normalized.
- `metric::`: Metric tensor at `position`.
- `spacetime::`: Spacetime.
- `model::`: Radiative model.
- `coords_top::`: Coordinates topology.
- `cache::`: Model cache.

# See also

- [`allocate_cache(model::AbstractRadiativeModel)`](@ref)
"""
function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top,
    cache)
    error("line_emission_profile not defined for this model.")
end
function is_final_position_at_source(position, spacetime, model)
    error("is_final_position_at_source not defined for this model.")
end

#Optional for surface emission models
"""
    allocate_cache(model::AbstractRadiativeModel)

Allocate a cache object for the given model. The cache object is used to store temporary data in radiative-model-related calculations.
"""
allocate_cache(::AbstractRadiativeModel) = nothing

function surface_differential!(differential, position, model, coords_top)
    error("Surface differential not defined for this model.")
end

"""
    temperature(position::AbstractVector, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)

    Temperature of the model at given `position` in the `spacetime`.
"""
temperature(position, spacetime, model) = error("Temperature not defined for this model.")

include("radiativeprocesses/thermalemission.jl")
include("radiativeprocesses/bremsstrahlung.jl")
include("radiativeprocesses/synchrotron.jl")

include("accretiondisk/accretiondisk.jl")
include("circularhotspot.jl")
include("onionhotspots.jl")
include("staracrosswormhole.jl")
include("verticalscreen.jl")
include("lamppostcorona.jl")
include("iontorus.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")

function unit_surface_normal!(vector, position, metric, metric_inverse, model, coords_top)
    surface_differential!(vector, position, model, coords_top)
    vector .= raise_index(vector, metric_inverse)
    normalize_spacelike!(vector, metric)
    return nothing
end

function rest_frame_four_velocity!(v,
    position,
    metric,
    spacetime,
    model,
    coords_top,
    ::Nothing)
    rest_frame_four_velocity!(v, position, metric, spacetime, model, coords_top)
end

function rest_frame_four_velocity!(v,
    position,
    metric,
    spacetime,
    model,
    coords_top,
    ::Nothing,
    ::Nothing)
    rest_frame_four_velocity!(v, position, metric, spacetime, model, coords_top)
end

function rest_frame_four_velocity!(v,
    position,
    metric,
    spacetime,
    model,
    coords_top,
    spacetime_cache::AbstractSpacetimeCache,
    ::Nothing)
    rest_frame_four_velocity!(v, position, metric, spacetime, model, coords_top, spacetime_cache)
end

function rest_frame_four_velocity!(v,
    position,
    metric,
    spacetime,
    model,
    coords_top,
    ::Nothing,
    model_cache::AbstractModelCache)
    rest_frame_four_velocity!(v, position, metric, spacetime, model, coords_top, model_cache)
end

function rest_frame_bolometric_intensity(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top,
    ::Nothing)
    rest_frame_bolometric_intensity(position,
        momentum,
        rest_frame_four_velocity,
        metric,
        spacetime,
        model,
        coords_top)
end
function rest_frame_specific_intensity(position,
    momentum,
    energy,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model,
    coords_top,
    ::Nothing)
    rest_frame_specific_intensity(position,
        momentum,
        energy,
        rest_frame_four_velocity,
        metric,
        spacetime,
        model,
        coords_top)
end


"""
    lorentz_factors(positions, spacetime, model)

    Lorentz factors of the rest frame four velocities of `model` at given `positions`.
"""
function lorentz_factors(positions::AbstractMatrix, 
    spacetime::AbstractSpacetime, 
    model::AbstractRadiativeModel) 
    coords_top = coordinates_topology(spacetime)
    γ = zeros(length(positions))
    g = zeros(4,4)
    u = zeros(4)
    for (i, position) in enumerate(positions)
        metric!(g, position, spacetime)
        rest_frame_four_velocity!(u, position, g, spacetime, model, coords_top)
        γ[i] = u[1] 
    end
    return γ
end