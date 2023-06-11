@kwdispatch geometrized_to_CGS
@kwdispatch CGS_to_geometrized

@kwmethod geometrized_to_CGS(magnitude, dim; M1::Number) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=M1) 
@kwmethod CGS_to_geometrized(magnitude, dim; M1::Number) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=M1)

@kwmethod geometrized_to_CGS(magnitude, dim, configurations::AbstractConfigurations;) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=configurations.unit_mass_in_solar_masses)
@kwmethod CGS_to_geometrized(magnitude, dim, configurations::AbstractConfigurations;) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=configurations.unit_mass_in_solar_masses)

function geometrized_unit_magnitude_in_CGS(dim; M1)
    #Returns the value in CGS of the geometrized system unit of dimension dim
    #M1 is the chosen unit mass in solar masses 

    c = PhysicalConstants.c
    G = PhysicalConstants.G
    M_sun = PhysicalConstants.M_sun

    M1_cgs = M1*M_sun        #Geometrized unit mass in CGS 
    L1_cgs = G*M1_cgs/c^2    #Geometrized unit length in CGS
    T1_cgs = L1_cgs/c        #Geometrized unit time in CGS

    return L1_cgs^dim.L*M1_cgs^dim.M*T1_cgs^dim.T
end

pc_to_cm(d) = d*PhysicalConstants.pc
cm_to_pc(d) = d/PhysicalConstants.pc
eV_to_erg(E) = E*PhysicalConstants.eV
erg_to_eV(E) = E/PhysicalConstants.eV
keV_to_erg(E) = E*(1e3*PhysicalConstants.eV)
erg_to_keV(E) = E/(1e3*PhysicalConstants.eV)
per_eV_to_per_erg(E) = erg_to_eV(E)
per_erg_to_per_eV(E) = eV_to_erg(E)
per_keV_to_per_erg(E) = erg_to_keV(E)
per_erg_to_per_keV(E) = keV_to_erg(E)