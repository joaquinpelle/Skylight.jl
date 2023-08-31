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
allocate_cache(::AbstractRadiativeModel) = nothing
function surface_differential!(differential, position, model, coords_top)
    error("Surface differential not defined for this model.")
end

#For accretion disks
temperature(position, spacetime, model) = error("Temperature not defined for this model.")

include("radiativeprocesses/thermalemission.jl")
include("radiativeprocesses/bremsstrahlung.jl")
include("radiativeprocesses/synchrotron.jl")

include("syntheticpolarcap.jl")
include("onionhotspots.jl")
include("bogdanovpolarcap.jl")
include("accretiondisk.jl")
include("staracrosswormhole.jl")
include("verticalscreen.jl")
include("lamppostcorona.jl")
include("iontorus.jl")
include("dummyextendedregion.jl")
include("dummymodel.jl")
include("../spacetimes/superposedpn/bbhdisk.jl")
include("../spacetimes/superposedpn/bbhcorona.jl")

stationarity(::AbstractRadiativeModel) = IsNotStationary()

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
    for i in axes(positions, 2)
        metric!(g, position, spacetime)
        rest_frame_four_velocity!(u, position, g, spacetime, model, coords_top)
        γ[i] = u[1] 
    end
    return γ
end

"""
    lorentz_factors(positions, spacetime, model)

    Lorentz factors of the rest frame four velocities of `model` at given `positions`.
"""
function lorentz_factors(positions::AbstractMatrix, 
    spacetime::AbstractSpacetime, 
    model::AbstractRadiativeModel,
    spacetime_cache::AbstractSpacetimeCache) 
    coords_top = coordinates_topology(spacetime)
    γ = zeros(length(positions))
    g = zeros(4,4)
    u = zeros(4)
    for i in axes(positions, 2)
        metric!(g, position, spacetime, spacetime_cache)
        rest_frame_four_velocity!(u, position, g, spacetime, model, coords_top, spacetime_cache)
        γ[i] = u[1] 
    end
    return γ
end