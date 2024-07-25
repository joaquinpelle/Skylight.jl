@doc raw"""
    NovikovThorneDisk <: AbstractAccretionDisk

Novikov & Thorne geometrically thin, optically thick accretion disk model around a Schwarzschild/Kerr black hole.

# Fields
- `inner_radius::Float64`: The inner radius of the accretion disk. Must be larger than or equal to zero.
- `outer_radius::Float64`: The outer radius of the accretion disk. Must be larger than or equal to `inner_radius`.
- `M1::Float64`: The unitary mass in solar masses. Must be positive.
- `Mdot_to_MEdd::Float64`: The accretion rate in units of the Eddington accretion rate. Must be positive.
- `η::Float64`: The radiative efficiency of the disk, which must be in the range (0, 1].
- `rotation_sense::AbstractRotationSense`: The sense of rotation of the disk, which can be either `ProgradeRotation()` or `RetrogradeRotation()`. Default is `ProgradeRotation()`.

# Examples
```julia
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M=1.0, a=0.5)
risco = isco_radius(spacetime, ProgradeRotation())
disk = NovikovThorneDisk(inner_radius = risco, outer_radius=1000.0, M1=1e7, Mdot_to_MEdd=0.1, η=0.1)
```
"""
@with_kw struct NovikovThorneDisk{T} <: AbstractAccretionDisk
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

function temperature(position, 
    spacetime::Union{AbstractSchwarzschildSpacetime,AbstractKerrSpacetime}, 
    model::NovikovThorneDisk)
    c2 = PhysicalConstants.c2
    σ = PhysicalConstants.σ
    M1 = model.M1
    Mdot_to_MEdd = model.Mdot_to_MEdd
    η = model.η

    r = radius(position, spacetime)
    M = mass(spacetime)
    a = spin(spacetime)
    χ = a/M
    risco = isco_radius(spacetime, model.rotation_sense)

    rCGS = geometrized_to_CGS(r, Dimensions.length, M1 = M1)
    MCGS = geometrized_to_CGS(M, Dimensions.mass, M1 = M1)
    Mdot = Mdot_to_MEdd*Eddington_accretion_rate(MCGS, η)

    x = sqrt(r / M)
    x0 = sqrt(risco/M)
    x1 = 2 * cos((acos(χ) - π) / 3)
    x2 = 2 * cos((acos(χ) + π) / 3)
    x3 = -2 * cos(acos(χ) / 3)

    f = 1.5/(x^3 - 3*x + 2χ)*(x - x0 - (3/2) * χ * log(x/x0) - 
              (3 * (x1 - χ)^2 / (x1 * (x1 - x2) * (x1 - x3))) * log((x - x1)/(x0 - x1)) -
              (3 * (x2 - χ)^2 / (x2 * (x2 - x1) * (x2 - x3))) * log((x - x2)/(x0 - x2)) -
              (3 * (x3 - χ)^2 / (x3 * (x3 - x1) * (x3 - x2))) * log((x - x3)/(x0 - x3)))

    F = Mdot*c2/(4*π*rCGS^2)*f
    T = (F/σ)^0.25
    return T
end

axisymmetry(::NovikovThorneDisk) = IsAxisymmetric()