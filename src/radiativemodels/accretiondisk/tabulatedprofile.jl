# Accretion disk with tabulated emission profile
@with_kw struct AccretionDiskWithTabulatedProfile{T, S} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation()
    filename::String
    profile_interpolator::S = build_interpolator(filename)
end

function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model::AccretionDiskWithTabulatedProfile,
    coords_top,
    cache)
    r = radius(position, spacetime)
    return model.profile_interpolator(r)
end
