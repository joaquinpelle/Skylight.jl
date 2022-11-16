export DummyModel

struct DummyModel <: RadiativeModel end

coordinate_system_class(DummyExtendedRegion) = CartesianClass()

function set_model_four_velocity!(vector, position, metric, model::DummyModel, coord_system)
    vector .= ∂t()
end
