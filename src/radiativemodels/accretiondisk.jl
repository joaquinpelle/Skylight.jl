function set_surface_differential!(covector, position, ::AbstractAccretionDisk, ::CartesianTopology)
    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 0.0
    covector[4] = 1.0
end

function set_surface_differential!(covector, position, ::AbstractAccretionDisk, ::SphericalTopology)
    covector[1] = 0.0
    covector[2] = 0.0
    covector[3] = 1.0
    covector[4] = 0.0
end    

function set_emitter_four_velocity!(vector, position, metric, spacetime, model::AbstractAccretionDisk, coords_top)
    angular_speed = circular_geodesic_angular_speed(position, spacetime, model.rotation_sense)
    tangent_vector_zaxis_rotation!(vector, position, angular_speed, metric, coords_top)
end

function emitted_bolometric_intensity(position, momentum, emitter_four_velocity, metric, spacetime, model::AbstractAccretionDisk, coords_top)
    T = temperature(position, spacetime, model)
    return thermal_emission_bolometric_intensity(T)
end

function emitted_specific_intensity(position, momentum, energy, emitter_four_velocity, metric, spacetime, model::AbstractAccretionDisk, coords_top)
    T = temperature(position, spacetime, model)
    return thermal_emission_specific_intensity(energy, T)
end

## The following function is used to check if the ray is inside the accretion disk
function is_final_position_at_source(position, spacetime, model::AbstractAccretionDisk)
    r = radius(position, spacetime)
    return (r >= model.inner_radius) && (r <= model.outer_radius)
end

# Novikov-Thorne disk

@with_kw struct NovikovThorneDisk{T} <: AbstractAccretionDisk 
    inner_radius::Float64
    outer_radius::Float64
    rotation_sense::T = ProgradeRotation() 
end

#Dummy function to be defined later
function temperature(position, spacetime, ::NovikovThorneDisk)
    r = radius(position, spacetime)
    return 1.5e-8 * (1.0 - 2.0 / r)^(0.25)
end

# Accretion disk with tabulated temperature

@with_kw struct AccretionDiskWithTabulatedTemperature{T,S} <: AbstractAccretionDisk
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
