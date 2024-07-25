@doc raw"""
    ShakuraSunyaevDisk <: AbstractAccretionDisk

Shakura & Sunyaev geometrically thin, optically thick accretion disk model.

# Fields
- `inner_radius::Float64`: The inner of the accretion disk. Must be larger than or equal to zero.
- `outer_radius::Float64`: The outer radius of the accretion disk. Must be larger than or equal to `inner_radius`.
- `M1::Float64`: The unitary mass in solar masses. Must be positive.
- `Mdot_to_MEdd::Float64`: The accretion rate in units of the Eddington accretion rate. Must be positive.
- `η::Float64`: The radiative efficiency of the disk, which must be in the range (0, 1].
- `rotation_sense::AbstractRotationSense`: The sense of rotation of the disk, which can be either `ProgradeRotation()` or `RetrogradeRotation()`. Default is `ProgradeRotation()`.

# Examples
```julia
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.5)
risco = isco_radius(spacetime, ProgradeRotation())
disk = ShakuraSunyaevDisk(inner_radius = risco, outer_radius=1000.0, M1=1e7, Mdot_to_MEdd=0.1, η=0.1)
```
"""
@with_kw struct ShakuraSunyaevDisk{T} <: AbstractAccretionDisk
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
    @assert 0.0<η<=1.0 "η must be in the range (0,1]"
end

function temperature(position, spacetime::AbstractBlackHoleSpacetime, model::ShakuraSunyaevDisk)
    σ = PhysicalConstants.σ
    G = PhysicalConstants.G
    M1 = model.M1
    Mdot_to_MEdd = model.Mdot_to_MEdd
    η = model.η

    r = radius(position, spacetime)
    M = mass(spacetime)
    rin = isco_radius(spacetime, model.rotation_sense)
    rCGS = geometrized_to_CGS(r, Dimensions.length, M1 = M1)
    MCGS = geometrized_to_CGS(M, Dimensions.mass, M1 = M1)
    Mdot = Mdot_to_MEdd*Eddington_accretion_rate(MCGS, η)
    f = 1-(rin/r)^0.5
    T = (3Mdot/(8*π*σ)*G*MCGS/rCGS^3*f)^0.25
    return T
end

axisymmetry(::ShakuraSunyaevDisk) = IsAxisymmetric()