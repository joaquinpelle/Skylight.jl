export DummyModel

struct DummyModel <: RadiativeModel end

function set_emitter_four_velocity!(vector, position, metric, spacetime, model::DummyModel, coord_system)
    vector .= ∂t()
end
