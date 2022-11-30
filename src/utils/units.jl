module Dimensions

using Constants

export geometrized_to_CGS
export CGS_to_geometrized

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

geometrized_to_CGS(magnitude; dim, M1) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=M1) 
CGS_to_geometrized(magnitude; dim, M1) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=M1) 

function geometrized_unit_magnitude_in_CGS(dim; M1)

    #Returns the value in CGS of the geometrized system unit of dimension dim
    #M1 is the geometrized system unit mass in solar masses 

    c = Constants.c
    G = Constants.G
    M_sun = Constants.M_sun

    M1_cgs = M1*M_sun        #Geometrized unit mass in CGS 
    L1_cgs = G*M1_cgs/c^2    #Geometrized unit length in CGS
    T1_cgs = L1_CGS/c        #Geometrized unit time in CGS

    return L1_cgs^dim.L*M1_cgs^dim.M*T1_cgs^dim.T

end

end