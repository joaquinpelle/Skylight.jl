include("shakurasunyaev.jl")
include("novikovthorne.jl")
include("rar.jl")
include("flatlamppostprofile.jl")
include("tabulatedprofile.jl")
include("tabulatedtemperature.jl")

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

#TODO note that it's assumed to be geom thin, opt thick, stationary
stationarity(::AbstractAccretionDisk) = IsStationary()

## The following function is used to check if the ray is inside the accretion disk
function is_final_position_at_source(position, spacetime, model::AbstractAccretionDisk)
    r = radius(position, spacetime)
    return (r >= model.inner_radius) && (r <= model.outer_radius)
end

radial_bins(disk::AbstractAccretionDisk; nbins) = range(disk.inner_radius, disk.outer_radius, length=nbins+1)

function Eddington_luminosity(M)
    return  PhysicalConstants.LEdd_sun*(M/PhysicalConstants.M_sun) 
end

function Eddington_accretion_rate(M, η)
    return Eddington_luminosity(M)/(η*PhysicalConstants.c^2)
end