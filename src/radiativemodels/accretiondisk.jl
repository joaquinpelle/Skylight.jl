function surface_differential!(covector,
    position,
    ::AbstractAccretionDisk,
    ::CartesianTopology)
    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 0.0
    covector[4] = 1.0
end

function surface_differential!(covector,
    position,
    ::AbstractAccretionDisk,
    ::SphericalTopology)
    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0
end

function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    model::AbstractAccretionDisk,
    coords_top)
    angular_speed = circular_geodesic_angular_speed(position,
        spacetime,
        model.rotation_sense)
    circular_motion_four_velocity!(vector, position, angular_speed, metric, coords_top)
end

function rest_frame_bolometric_intensity(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model::AbstractAccretionDisk,
    coords_top)
    T = temperature(position, spacetime, model)
    return thermal_emission_bolometric_intensity(T)
end

function rest_frame_specific_intensity(position,
    momentum,
    energy,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model::AbstractAccretionDisk,
    coords_top)
    T = temperature(position, spacetime, model)
    return thermal_emission_specific_intensity(energy, T)
end

stationarity(::AbstractAccretionDisk) = IsStationary()

## The following function is used to check if the ray is inside the accretion disk
function is_final_position_at_source(position, spacetime, model::AbstractAccretionDisk)
    r = radius(position, spacetime)
    return (r >= model.inner_radius) && (r <= model.outer_radius)
end

radial_bins(disk::AbstractAccretionDisk; nbins) = range(disk.inner_radius, disk.outer_radius, length=nbins+1)

# Novikov-Thorne disk
@with_kw struct NovikovThorneDisk{T} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation()

    @assert inner_radius>=0.0 "Inner radius must be non-negative"
    @assert outer_radius>=inner_radius "Outer radius must be larger than inner radius"
    @assert isa(rotation_sense, AbstractRotationSense) "Rotation sense must be either ProgradeRotation() or RetrogradeRotation()"
end

#Dummy function to be defined later
function temperature(position, spacetime, ::NovikovThorneDisk)
    r = radius(position, spacetime)
    return 1.5e-8 * (1.0 - 2.0 / r)^(0.25)
end

function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime::KerrSpacetimeBoyerLindquistCoordinates,
    ::NovikovThorneDisk,
    coords_top,
    cache)
    r = radius(position, spacetime)
    return 1.0 / r^2
end

# Shakura-Sunyaev disk

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

# RAR disk
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

# Accretion disk with tabulated temperature
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

# Accretion disk with tabulated emission profile
@with_kw struct AccretionDiskWithTabulatedProfile{T, S} <: AbstractAccretionDisk
    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation()
    filename::String
    profile_interpolator::S = build_interpolator(filename)
end

function line_emission_profile(position,
    momentum,
    rest_frame_four_velocity,
    metric,
    spacetime,
    model::AccretionDiskWithTabulatedProfile,
    coords_top,
    cache)
    r = radius(position, spacetime)
    return model.profile_interpolator(r)
end

# Accretion disk with tabulated emission profile
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

function Eddington_luminosity(M)
    return 4π*PhysicalConstants.c*PhysicalConstants.G*M*PhysicalConstants.mp/PhysicalConstants.sigma_T
end

function Eddington_accretion_rate(M, η)
    return Eddington_luminosity(M)/(η*PhysicalConstants.c^2)
end