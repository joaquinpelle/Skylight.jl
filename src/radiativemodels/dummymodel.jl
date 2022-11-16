export DummyModel

struct DummyModel <: RadiativeModel end

function set_model_four_velocity!(vector, position, metric, model::DummyModel, coord_system)
    vector .= ∂t()
end
