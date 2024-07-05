@doc raw"""
    AccretionDiskWithTabulatedProfile <: AbstractAccretionDisk

Geometrically thin, optically thick accretion disk model with tabulated line emission radial profile. This is useful, for example, for incorporating
the emissivity profiles calculated from a corona model.

# Fields
- `inner_radius::Float64`: The inner of the accretion disk. Must be larger than or equal to zero.
- `outer_radius::Float64`: The outer radius of the accretion disk. Must be larger than or equal to `inner_radius`.
- `rotation_sense::AbstractRotationSense`: The sense of rotation of the disk, which can be either `ProgradeRotation()` or `RetrogradeRotation()`. Default is `ProgradeRotation()`.
- `filename::String`: The name of a two-column file containing the line emission radial profile as radius vs. profile.

# Examples
```julia
disk = AccretionDiskWithTabulatedProfile(inner_radius=6.0, outer_radius=1000.0, filename="emissivity.txt")
```
"""
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
