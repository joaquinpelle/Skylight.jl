struct DummyModel <: AbstractRadiativeModel end

function set_emitter_four_velocity!(vector, position, metric, spacetime, model::DummyModel, coords_top)
    vector .= ∂t()
    return nothing
end
