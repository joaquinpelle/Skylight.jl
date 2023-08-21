module Dimensions

using Parameters

@with_kw struct Dimension{S}
    L::S = 0
    M::S = 0
    T::S = 0
end

function Base.:*(dim1::Dimension, dim2::Dimension)
    Dimension(L = dim1.L + dim2.L,
        M = dim1.M + dim2.M,
        T = dim1.T + dim2.T)
end

function Base.:/(dim1::Dimension, dim2::Dimension)
    Dimension(L = dim1.L - dim2.L,
        M = dim1.M - dim2.M,
        T = dim1.T - dim2.T)
end

Base.inv(dim::Dimension) = Dimension(L = -dim.L,
    M = -dim.M,
    T = -dim.T)

Base.:^(dim::Dimension, n::Number) = Dimension(L = dim.L * n,
    M = dim.M * n,
    T = dim.T * n)

const adim = Dimension()
const length = Dimension(L = 1)
const mass = Dimension(M = 1)
const time = Dimension(T = 1)

const action = length^2 * mass / time
const angle = adim
const angular_acceleration = time^-2
const angular_momentum = length^2 * mass / time
const angular_velocity = angle / time
const area = length^2
const area_density = length^-2 * mass
const energy = length^2 * mass / time^2
const energy_density = length^-1 * mass / time^2
const force = length * mass / time^2
const frequency = time^-1
const heat = length^2 * mass / time^2
const impulse = length * mass / time
const linear_density = length^-1 * mass
const mass_density = length^-3 * mass
const moment_of_inertia = length^2 * mass
const momentum = length * mass / time
const power = length^2 * mass / time^3
const pressure = length^-1 * mass / time^2
const specific_intensity = mass / time^2
const torque = length^2 * mass / time^2
const velocity = length / time
const volume = length^3
const wavenumber = length^-1

#Electromagnetical (Gaussian or Lorentz-Heaviside)
const capacitance = length
const charge = length^1.5 * mass^0.5 / time
const charge_density = length^-1.5 * mass^0.5 / time
const conductivity = time^-1
const current_density = length^-0.5 * mass^0.5 * time^-2
const current_intensity = length^1.5 * mass^0.5 * time^-2
const electric_potential = length^0.5 * mass^0.5 / time
const electric_field = length^-0.5 * mass^0.5 / time
const electric_displacement_field = length^-0.5 * mass^0.5 / time
const inductance = length^-1 * time^2
const magnetic_field = length^-0.5 * mass^0.5 / time
const magnetic_H_field = length^-0.5 * mass^0.5 / time
const magnetic_diffusivity = length^2 / time
const magnetic_dipole_moment = length^2.5 * mass^0.5 / time
const magnetic_flux = length^1.5 * mass^0.5 / time
const magnetization = length^0.5 * mass^0.5 / time^2
const permeability = length^-2 * time^2
const permittivity = adim
const resistance = length^-1 * time
const resistivity = time

end
