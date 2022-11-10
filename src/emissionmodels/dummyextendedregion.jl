struct DummyExtendedRegion <: EmissionModel end

coordinate_system_class(DummyExtendedRegion) = CartesianClass()

function set_model_four_velocity!(vector, position, metric, model::DummyExtendedRegion, coord_system)
    vector .= ∂t()
end