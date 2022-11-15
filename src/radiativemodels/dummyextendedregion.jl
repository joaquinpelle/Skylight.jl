export DummyExtendedRegion

@with_kw struct DummyExtendedRegion <: RadiativeModel 

    rbound::Float64 = 0.2

end

coordinate_system_class(DummyExtendedRegion) = SphericalClass()

function set_model_four_velocity!(vector, position, metric, model::DummyExtendedRegion, coord_system)
    vector .= ∂t()
end

function set_invariant_absorptivity!(αε, point, ε, model::DummyExtendedRegion)

    t, r, θ, φ = point
    @. αε = ε*exp(-(r-4)^2/2)*exp(-ε/2)

end

function set_invariant_emissivity!(jε, point, ε, model::DummyExtendedRegion)

    t, r, θ, φ = point
    @. jε = exp(-(r-6)^2/2)*exp(-ε/2)/ε^2

end

