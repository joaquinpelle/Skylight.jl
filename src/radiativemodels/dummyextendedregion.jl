struct DummyExtendedRegion <: AbstractRadiativeModel end

function emitter_four_velocity!(vector, position, metric, spacetime, model::DummyExtendedRegion, coords_top)
    vector .= ∂t()
    return nothing
end

function invariant_absorptivity!(αε, point, ε, model::DummyExtendedRegion)

    t, r, θ, φ = point
    @. αε = ε*exp(-(r-4)^2/2)*exp(-ε/2)

    return nothing
end

function invariant_emissivity!(jε, point, ε, model::DummyExtendedRegion)

    t, r, θ, φ = point
    @. jε = exp(-(r-6)^2/2)*exp(-ε/2)/ε^2

    return nothing
end

