@with_kw struct AccretionDiskWithFlatLamppostProfile{T} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    corona_height::Float64
    rotation_sense::T = ProgradeRotation()
end

function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model::AccretionDiskWithFlatLamppostProfile,
    coords_top,
    cache)
    r = radius(position, spacetime)
    h = model.corona_height
    return h/(r^2+h^2)^1.5
end
