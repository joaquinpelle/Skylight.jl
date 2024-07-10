@doc raw"""
    rest_frame_four_velocity!(vector::AbstractVector, position::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, spacetime_cache::AbstractSpacetimeCache, model_cache::AbstractModelCache)

    Rest frame four velocity of the model at given `position`. This is the frame where the model radiative functions are defined.

# Arguments

- `vector`: Output vector.
- `position`: Position where the four-velocity is evaluated.
- `metric`: Metric tensor at `position`.
- `spacetime`: Spacetime.
- `model`: Radiative model.
- `coords_top`: Coordinates topology.
- `spacetime_cache`: Spacetime cache.
- `model_cache`: Model cache.

# See also

- [`allocate_cache(model::AbstractRadiativeModel)`](@ref)
"""
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

@doc raw"""
    rest_frame_bolometric_intensity(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)

    Bolometric intensity of the model at given `position`.

# Arguments

- `position`: Position where the intensity is evaluated.
- `momentum`: Momentum of the emission (frequency and direction).
- `rest_frame_four_velocity`: Rest frame four velocity of the model at `position`. Must but be normalized.
- `metric`: Metric tensor at `position`.
- `spacetime`: Spacetime.
- `model`: Radiative model.
- `coords_top`: Coordinates topology.
- `cache`: Model cache.

# See also

- [`allocate_cache(model::AbstractRadiativeModel)`](@ref)
"""
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

@doc raw"""
    rest_frame_specific_intensity(position::AbstractVector, momentum::AbstractVector, energy::Real, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)

    Specific intensity of the model at given `position`.

# Arguments

- `position`: Position where the intensity is evaluated.
- `momentum`: Momentum of the emission (frequency and direction).
- `energy`: Energy of the emission. This should be assumed to be in CGS units.
- `rest_frame_four_velocity`: Rest frame four velocity of the model at `position`. Must but be normalized.
- `metric`: Metric tensor at `position`.
- `spacetime`: Spacetime.
- `model`: Radiative model.
- `coords_top`: Coordinates topology.
- `cache`: Model cache.

# See also

- [`allocate_cache(model::AbstractRadiativeModel)`](@ref)
"""
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

@doc raw"""
    line_emission_profile(position::AbstractVector, momentum::AbstractVector, rest_frame_four_velocity::AbstractVector, metric::AbstractMatrix, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)

    Emissivity profile for line emission radiative models.

# Arguments

- `position`: Position of the emission.
- `momentum`: Momentum of the emission (frequency and direction). Must be null.
- `rest_frame_four_velocity`: Rest frame four velocity of the model at `position`. Must but be normalized.
- `metric`: Metric tensor at `position`.
- `spacetime`: Spacetime.
- `model`: Radiative model.
- `coords_top`: Coordinates topology.
- `cache`: Model cache.

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

@doc raw"""
    is_final_position_at_source(position::AbstractVector, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)

    Check if the final position of the photon is at the source of the model. This function is used in the observer-to-emitter method
    in vacuum to discard rays that do not intersect the source.

# Arguments

- `position`: Final position of the geodesic.
- `spacetime`: Spacetime.
- `model`: Radiative model.
"""
function is_final_position_at_source(position, spacetime, model)
    error("is_final_position_at_source not defined for this model.")
end

#Optional for surface emission models
@doc raw"""
    allocate_cache(model::AbstractRadiativeModel)

Allocate a cache object for the given model. The cache object is used to store temporary data in radiative-model-related calculations.
"""
allocate_cache(::AbstractRadiativeModel) = nothing

@doc raw"""
    surface_differential!(differential::AbstractVector, position::AbstractVector, model::AbstractSurfaceEmissionModel, coords_top::AbstractCoordinatesTopology)

    Differential of the function defining the emitting surface in the model. For example, for an emitting sphere in Cartesian coordinates, the output would be `[0,2x,2y,2z]`.
    The normalization of the differential is not important, as it is only used to calculate the surface (unit) normal in terms of the metric. 

# Arguments

- `differential`: Output vector.
- `position`: Position where the differential is evaluated.
- `model`: Radiative model.
- `coords_top`: Coordinates topology.
"""
function surface_differential!(differential, position, model, coords_top)
    error("Surface differential not defined for this model.")
end

@doc raw"""
    temperature(position::AbstractVector, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)

    Temperature of the model at given `position`.
"""
temperature(position, spacetime, model) = error("Temperature not defined for this model.")

@doc raw"""
    space_positions(npoints::Integer, spacetime::AbstractSpacetime, model::AbstractRadiativeModel, coords_top::AbstractCoordinatesTopology, cache::AbstractModelCache)

    Generate `npoints` space (3D) positions where to launch photon packages form in the EtO method. The positions are generated in the model's coordinates and should be used to calculate the model's specific intensity.


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


@doc raw"""
    lorentz_factors(positions, spacetime::AbstractSpacetime, model::AbstractRadiativeModel)

    Lorentz factors of the rest frame four velocities of `model` at given list of positions. Positions should be an iterable
    object with the positions as elements.

    # Returns

    An array with the Lorentz factors of the rest frame four velocities of the model at the given positions.

    # Example

    ```
    positions = [[0.0, 5.0, 0.0, 0.0], [1.0, 0.0, 5.0, 0.0]]
    spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
    model = DummyModel()
    γ = lorentz_factors(positions, spacetime, model)
    ```
"""
function lorentz_factors(positions, 
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