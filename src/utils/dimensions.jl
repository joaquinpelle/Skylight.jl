module Dimensions

using Parameters

@with_kw struct Dimension{S}
    L::S = 0
    M::S = 0
    T::S = 0
end

Base.:*(dim1::Dimension, dim2::Dimension) = Dimension(L=dim1.L+dim2.L, 
                                                      M=dim1.M+dim2.M, 
                                                      T=dim1.T+dim2.T)

Base.:/(dim1::Dimension, dim2::Dimension) = Dimension(L=dim1.L-dim2.L, 
                                                      M=dim1.M-dim2.M, 
                                                      T=dim1.T-dim2.T)

Base.inv(dim::Dimension) = Dimension(L=-dim.L, 
                                     M=-dim.M, 
                                     T=-dim.T)

Base.:^(dim::Dimension, n::Number) = Dimension(L=dim.L*n,
                                                M=dim.M*n,
                                                T=dim.T*n)

adim = Dimension()                                                
length = Dimension(L=1)
mass = Dimension(M=1)
time = Dimension(T=1)

action = length^2*mass/time
angle = adim
angular_acceleration = time^-2
angular_momentum = length^2*mass/time
angular_velocity = angle/time
area = length^2
area_density = length^-2*mass
energy = length^2*mass/time^2
energy_density = length^-1*mass/time^2
force = length*mass/time^2
frequency = time^-1
heat = length^2*mass/time^2
impulse = length*mass/time
linear_density = length^-1*mass
mass_density = length^-3*mass
moment_of_inertia = length^2*mass
momentum = length*mass/time
power = length^2*mass/time^3
pressure = length^-1*mass/time^2
specific_intensity = mass/time^2
torque = length^2*mass/time^2
velocity = length/time
volume = length^3
wavenumber = length^-1

#Electromagnetical (Gaussian or Lorentz-Heaviside)
capacitance = length
charge = length^1.5*mass^0.5/time
charge_density = length^-1.5*mass^0.5/time
conductivity = time^-1
current_density = length^-0.5*mass^0.5*time^-2
current_intensity = length^1.5*mass^0.5*time^-2
electric_potential = length^0.5*mass^0.5/time
electric_field = length^-0.5*mass^0.5/time
electric_displacement_field = length^-0.5*mass^0.5/time
inductance = length^-1*time^2
magnetic_field = length^-0.5*mass^0.5/time
magnetic_H_field = length^-0.5*mass^0.5/time
magnetic_diffusivity = length^2/time
magnetic_dipole_moment = length^2.5*mass^0.5/time
magnetic_flux = length^1.5*mass^0.5/time
magnetization = length^0.5*mass^0.5/time^2
permeability = length^-2*time^2
permittivity = adim
resistance = length^-1*time
resistivity = time

end