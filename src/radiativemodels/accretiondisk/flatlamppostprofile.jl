@doc raw"""
    AccretionDiskWithFlatLamppostProfile <: AbstractAccretionDisk

Geometrically thin, optically thick accretion disk model with lamppost corona line emissivity profile as in flat spacetime:

``\epsilon \propto \frac{h}{(r^2+h^2)^{3/2}}``

where `h` is the height of the corona above the disk plane.

# Fields
- `inner_radius::Float64`: The inner of the accretion disk. Must be larger than or equal to zero.
- `outer_radius::Float64`: The outer radius of the accretion disk. Must be larger than or equal to `inner_radius`.
- `rotation_sense::AbstractRotationSense`: The sense of rotation of the disk, which can be either `ProgradeRotation()` or `RetrogradeRotation()`. Default is `ProgradeRotation()`.
- `corona_height::String`: The height of the corona above the disk plane.

# Examples
```julia
disk = AccretionDiskWithFlatLamppostProfile(inner_radius=6.0, outer_radius=1000.0, corona_height=5.0)
```
"""
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
