@with_kw struct DexterAgolDisk{T} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation()

    @assert inner_radius>=0.0 "Inner radius must be non-negative"
    @assert outer_radius>=inner_radius "Outer radius must be larger than inner radius"
    @assert isa(rotation_sense, AbstractRotationSense) "Rotation sense must be either ProgradeRotation() or RetrogradeRotation()"
end

function temperature(position, spacetime, ::DexterAgolDisk)
    r = radius(position, spacetime)
    return 1.5e-8 * (1.0 - 2.0 / r)^(0.25)
end

function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime::AbstractKerrSpacetime,
    model::DexterAgolDisk,
    coords_top,
    cache)
    r = radius(position, spacetime)
    return 1.0/r^2
end

axisymmetry(::DexterAgolDisk) = IsAxisymmetric()