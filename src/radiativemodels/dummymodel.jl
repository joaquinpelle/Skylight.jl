struct DummyModel <: AbstractRadiativeModel end

stationarity(::DummyModel) = IsStationary()
spherical_symmetry(::DummyModel) = IsSphericallySymmetric()

function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    model::DummyModel,
    coords_top)
    vector .= ∂t()
    return nothing
end
