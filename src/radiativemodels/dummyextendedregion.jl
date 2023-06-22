struct DummyExtendedRegion <: AbstractRadiativeModel end

function emitter_four_velocity!(vector, position, metric, spacetime, model::DummyExtendedRegion, coords_top)
    vector .= ∂t()
    return nothing
end

function invariant_absorptivity!(αε, position, ε, ::DummyExtendedRegion)

    r = position[2]
    @. αε = ε*exp(-(r-4)^2/2)*exp(-ε/2)

    return nothing
end

function invariant_emissivity!(jε, position, ε, ::DummyExtendedRegion)

    r = position[2]
    @. jε = exp(-(r-6)^2/2)*exp(-ε/2)/ε^2

    return nothing
end

