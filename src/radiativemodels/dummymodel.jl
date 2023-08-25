struct DummyModel <: AbstractRadiativeModel end

function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    model::DummyModel,
    coords_top)
    vector .= ∂t()
    return nothing
end
