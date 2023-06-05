struct DummyModel <: AbstractRadiativeModel end

function emitter_four_velocity!(vector, position, metric, spacetime, model::DummyModel, coords_top)
    vector .= ∂t()
    return nothing
end
