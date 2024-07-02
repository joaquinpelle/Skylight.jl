@doc raw"""
    AccretionDiskWithTabulatedTemperature <: AbstractAccretionDisk

Geometrically thin, optically thick accretion disk model with tabulated temperature.

# Fields
- `inner_radius::Float64`: The inner of the accretion disk. Must be larger than or equal to zero.
- `outer_radius::Float64`: The outer radius of the accretion disk. Must be larger than or equal to `inner_radius`.
- `rotation_sense::AbstractRotationSense`: The sense of rotation of the disk, which can be either `ProgradeRotation()` or `RetrogradeRotation()`. Default is `ProgradeRotation()`.
- `filename::String`: The name of a two-column file containing the temperature profile as radius vs. temperature.

# Examples
```julia
disk = AccretionDiskWithTabulatedTemperature(inner_radius=6.0, outer_radius=1000.0, filename="temperature.txt")
```
"""
@with_kw struct AccretionDiskWithTabulatedTemperature{T, S} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation()
    filename::String
    temperature_interpolator::S = build_interpolator(filename)
end

function temperature(position, spacetime, model::AccretionDiskWithTabulatedTemperature)
    r = radius(position, spacetime)
    return model.temperature_interpolator(r)
end