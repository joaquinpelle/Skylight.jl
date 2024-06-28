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

function temperature(position, spacetime::AbstractKerrSpacetime, model::NovikovThorneDisk)
    M1 = model.M1
    Mdot_to_MEdd = model.Mdot_to_MEdd
    η = model.η
    c2 = PhysicalConstants.c2

    r = radius(position, spacetime)
    M = mass(spacetime)
    a = spin(spacetime)
    χ = a / M
    risco = innermost_circular_orbit_radius(spacetime, model.rotation_sense)

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