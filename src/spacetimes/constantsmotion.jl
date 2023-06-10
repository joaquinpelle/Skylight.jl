time_translation_killing_vector(spacetime::AbstractSpacetime) = _time_translation_killing_vector(stationarity(spacetime))
axial_rotation_killing_vector(position, spacetime::AbstractSpacetime) = _axial_rotation_killing_vector(position, coordinates_topology(spacetime), axial_symmetry(spacetime))
spherical_rotation_killing_vectors(position, spacetime::AbstractSpacetime) = _spherical_rotation_killing_vectors(position, coordinates_topology(spacetime), spherical_symmetry(spacetime))

_time_translation_killing_vector(::IsNotStationary) = error("Spacetime is not stationary. Time translation Killing vector not defined.")
_axial_rotation_killing_vector(position, coords_top, ::IsNotAxiallySymmetric) = error("Spacetime is not axially symmetric. Axial rotation Killing vector not defined.")
_spherical_rotation_killing_vectors(position, coords_top, ::IsNotSphericallySymmetric) = error("Spacetime is not spherically symmetric. Spherical rotation Killing vectors not defined.")

_time_translation_killing_vector(::IsStationary) = time_translation_generator()
_axial_rotation_killing_vector(position, coords_top, ::IsAxiallySymmetric) = axial_rotation_generator(position,coords_top) 
_spherical_rotation_killing_vectors(position, coords_top, ::IsSphericallySymmetric) = spherical_rotation_generators(position, coords_top)

function energy(position::AbstractVector, momentum::AbstractVector, spacetime::AbstractSpacetime)
    g = metric(position, spacetime)
    k = time_translation_killing_vector(spacetime)
    return -scalar_product(k, momentum, g)
end

function axial_angular_momentum(position::AbstractVector, momentum::AbstractVector, spacetime::AbstractSpacetime)
    g = metric(position, spacetime)
    k = axial_rotation_killing_vector(position, spacetime)
    return -scalar_product(k, momentum, g)
end

function angular_momentum(position::AbstractVector, momentum::AbstractVector, spacetime::AbstractSpacetime)
    g = metric(position, spacetime)
    k = spherical_rotation_killing_vectors(position, spacetime)
    return [-scalar_product(k[:,i], momentum, g) for i in 1:3]
end

function energy(position::AbstractVector, momenta::AbstractMatrix, spacetime::AbstractSpacetime)
    g = metric(position, spacetime)
    k = time_translation_killing_vector(spacetime)
    return [-scalar_product(k, momenta[:,i], g) for i in axes(momenta,2)]
end

function axial_angular_momentum(position::AbstractVector, momenta::AbstractMatrix, spacetime::AbstractSpacetime)
    g = metric(position, spacetime)
    k = axial_rotation_killing_vector(position, spacetime)
    return [-scalar_product(k, momenta[:,i], g) for i in axes(momenta,2)]
end

function angular_momentum(position::AbstractVector, momenta::AbstractMatrix, spacetime::AbstractSpacetime)
    g = metric(position, spacetime)
    k = spherical_rotation_killing_vectors(position, spacetime)
    return [[-scalar_product(k[:,j], momenta[:,i], g) for j in 1:3] for i in axes(momenta,2)]
end

function energy(positions::AbstractMatrix, momenta::AbstractMatrix, spacetime::AbstractSpacetime)
    @assert size(positions) == size(momenta) "Positions and momenta must have the same size."
    return [energy(positions[:,i], momenta[:,i], spacetime) for i in axes(momenta,2)]
end

function axial_angular_momentum(positions::AbstractMatrix, momenta::AbstractMatrix, spacetime::AbstractSpacetime)
    @assert size(positions) == size(momenta) "Positions and momenta must have the same size."
    return [axial_angular_momentum(positions[:,i], momenta[:,i], spacetime) for i in axes(momenta,2)]
end

function angular_momentum(positions::AbstractMatrix, momenta::AbstractMatrix, spacetime::AbstractSpacetime)
    @assert size(positions) == size(momenta) "Positions and momenta must have the same size."
    return [angular_momentum(positions[:,i], momenta[:,i], spacetime) for i in axes(momenta,2)]
end