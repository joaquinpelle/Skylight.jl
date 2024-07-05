@doc raw"""
    RARDisk <: AbstractAccretionDisk

RAR extension of the Shakura & Sunyaev geometrically thin, optically thick accretion disk model ([Millauro et al. 2024](https://www.aanda.org/articles/aa/abs/2024/05/aa48461-23/aa48461-23.html)).

# Fields
- `inner_radius::Float64`: The inner of the accretion disk. Must be larger than or equal to zero.
- `outer_radius::Float64`: The outer radius of the accretion disk. Must be larger than or equal to `inner_radius`.
- `M1::Float64`: The unitary mass in solar masses. Must be positive.
- `Mdot_to_MEdd::Float64`: The accretion rate in units of the Eddington accretion rate. Must be positive.
- `η::Float64`: The radiative efficiency of the disk, which must be in the range (0, 1].
- `rotation_sense::AbstractRotationSense`: The sense of rotation of the disk, which can be either `ProgradeRotation()` or `RetrogradeRotation()`. Default is `ProgradeRotation()`.

# Examples
```julia
disk = RARDisk(inner_radius=0.0, outer_radius=1000.0, M1=1e7, Mdot_to_MEdd=0.1, η=0.1)
```
"""
@with_kw struct RARDisk{T} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    M1::Float64
    Mdot_to_MEdd::Float64
    η::Float64
    rotation_sense::T = ProgradeRotation()

    @assert inner_radius>=0.0 "Inner radius must be non-negative"
    @assert outer_radius>=inner_radius "Outer radius must be larger than inner radius"
    @assert isa(rotation_sense, AbstractRotationSense) "Rotation sense must be either ProgradeRotation() or RetrogradeRotation()"
    @assert M1>0.0 "M1 must be positive"
    @assert Mdot_to_MEdd>0.0 "Mdot_to_MEdd must be positive"
    @assert η>0.0 "η must be positive"
    @assert η<=1.0 "η must be less than or equal to 1"
end

function temperature(position, spacetime, model::RARDisk)
    σ = PhysicalConstants.σ
    G = PhysicalConstants.G
    rin = model.inner_radius
    M1 = model.M1
    Mdot_to_MEdd = model.Mdot_to_MEdd
    η = model.η

    r = radius(position, spacetime)
    M = mass_enclosed(r, spacetime)
    Min = mass_enclosed(rin, spacetime)
    dM = mass_enclosed_derivative(r, spacetime)

    rCGS = geometrized_to_CGS(r, Dimensions.length, M1 = M1)
    MCGS = geometrized_to_CGS(M, Dimensions.mass, M1 = M1)
    Mdot = Mdot_to_MEdd*Eddington_accretion_rate(MCGS, η)
    f = 1-(Min*rin/(M*r))^0.5
    g = (1-r/(3M)*dM)
    T = (3Mdot/(8*π*σ)*G*MCGS/rCGS^3*f*g)^0.25
    return T
end

function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    ::RARDisk,
    coords_top,
    cache)
    r = radius(position, spacetime)
    return 1.0 / r^2
end

axisymmetry(::RARDisk) = IsAxisymmetric()