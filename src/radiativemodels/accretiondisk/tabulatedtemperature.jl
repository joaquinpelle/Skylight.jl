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
