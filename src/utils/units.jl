@kwdispatch geometrized_to_CGS
@kwdispatch CGS_to_geometrized

@kwmethod function geometrized_to_CGS(magnitude, dim; M1::Number)
    magnitude * geometrized_unitary_magnitude_in_CGS(dim; M1 = M1)
end
@kwmethod function CGS_to_geometrized(magnitude, dim; M1::Number)
    magnitude / geometrized_unitary_magnitude_in_CGS(dim; M1 = M1)
end

@kwmethod function geometrized_to_CGS(magnitude,
    dim,
    configurations::AbstractConfigurations;)
    magnitude *
    geometrized_unitary_magnitude_in_CGS(dim; M1 = configurations.unit_mass_in_solar_masses)
end
@kwmethod function CGS_to_geometrized(magnitude,
    dim,
    configurations::AbstractConfigurations;)
    magnitude /
    geometrized_unitary_magnitude_in_CGS(dim; M1 = configurations.unit_mass_in_solar_masses)
end

function geometrized_unitary_magnitude_in_CGS(dim; M1)
    #Returns the value in CGS of the geometrized system unit of dimension dim
    #M1 is the chosen unit mass in solar masses 

    c = PhysicalConstants.c
    G = PhysicalConstants.G
    M_sun = PhysicalConstants.M_sun

    M1_cgs = M1 * M_sun        #Geometrized unit mass in CGS 
    L1_cgs = G * M1_cgs / c^2    #Geometrized unit length in CGS
    T1_cgs = L1_cgs / c        #Geometrized unit time in CGS

    return L1_cgs^dim.L * M1_cgs^dim.M * T1_cgs^dim.T
end

pc_to_cm(d) = d * PhysicalConstants.pc
cm_to_pc(d) = d / PhysicalConstants.pc
kpc_to_cm(d) = d * (1e3*PhysicalConstants.pc)
cm_to_kpc(d) = d / (1e3*PhysicalConstants.pc)
Hz_to_erg(ν) = ν * PhysicalConstants.h
erg_to_Hz(E) = E / PhysicalConstants.h
eV_to_erg(E) = E * PhysicalConstants.eV
erg_to_eV(E) = E / PhysicalConstants.eV
keV_to_erg(E) = E * (1e3 * PhysicalConstants.eV)
erg_to_keV(E) = E / (1e3 * PhysicalConstants.eV)
per_eV_to_per_erg(E) = erg_to_eV(E)
per_erg_to_per_eV(E) = eV_to_erg(E)
per_keV_to_per_erg(E) = erg_to_keV(E)
per_erg_to_per_keV(E) = keV_to_erg(E)
eV_to_K(T) = eV_to_erg(T) / PhysicalConstants.k_B
K_to_eV(T) = erg_to_eV(T) * PhysicalConstants.k_B 
keV_to_K(T) = keV_to_erg(T) / PhysicalConstants.k_B
K_to_keV(T) = erg_to_keV(T) * PhysicalConstants.k_B 
deg_to_as(θ) = θ * 3600.0
as_to_deg(θ) = θ / 3600.0
deg_to_mas(θ) = θ * 3600.0*1e3
mas_to_deg(θ) = θ / (3600.0*1e3)