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
    model::RARDisk,
    coords_top,
    cache)
    r = radius(position, spacetime)
    return 1.0 / r^2
end