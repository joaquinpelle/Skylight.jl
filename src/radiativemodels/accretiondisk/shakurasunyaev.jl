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
    @assert η>0.0 "η must be positive"
    @assert η<=1.0 "η must be less than or equal to 1"
end

function temperature(position, spacetime, model::ShakuraSunyaevDisk)
    σ = PhysicalConstants.σ
    G = PhysicalConstants.G
    rin = model.inner_radius
    M1 = model.M1
    Mdot_to_MEdd = model.Mdot_to_MEdd
    η = model.η

    r = radius(position, spacetime)
    M = mass(spacetime)
    rCGS = geometrized_to_CGS(r, Dimensions.length, M1 = M1)
    MCGS = geometrized_to_CGS(M, Dimensions.mass, M1 = M1)
    Mdot = Mdot_to_MEdd*Eddington_accretion_rate(MCGS, η)
    f = 1-(rin/r)^0.5
    T = (3Mdot/(8*π*σ)*G*MCGS/rCGS^3*f)^0.25
    return T
end