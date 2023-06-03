energy(position, momentum, spacetime::AbstractSpacetime) = energy(position, momentum, spacetime, is_stationary(spacetime))  
angular_momentum(position, momentum, spacetime::AbstractSpacetime) = angular_momentum(position, momentum, spacetime, is_spherically_symmetric(spacetime))  
z_angular_momentum(position, momentum, spacetime::AbstractSpacetime) = z_angular_momentum(position, momentum, spacetime, is_axially_symmetric(spacetime))  

energy(position, momentum, spacetime,::IsNotStationary) = error("Spacetime is not stationary. Energy not defined.")
angular_momentum(position, momentum, spacetime, ::IsNotSphericallySymmetric) = error("Spacetime is not spherically symmetric. Angular momentum not defined. If the spacetime is axially symmetric you can use the z angular momentum")
z_angular_momentum(position, momentum, spacetime, ::IsNotAxiallySymmetric) = error("Spacetime is not axially symmetric. z angular momentum not defined.")

