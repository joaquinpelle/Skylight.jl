struct DummyExtendedRegion <: EmissionModel end

function set_model_four_velocity!(vector, position, metric, model::DummyExtendedRegion, coord_system)
    vector .= ∂t()
end

get_number_of_points(model::DummyExtendedRegion) = 3