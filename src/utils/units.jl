export geometrized_to_CGS
export CGS_to_geometrized

@kwdispatch geometrized_to_CGS
@kwdispatch CGS_to_geometrized

@kwmethod geometrized_to_CGS(magnitude, dim; M1::Number) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=M1) 
@kwmethod CGS_to_geometrized(magnitude, dim; M1::Number) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=M1)

@kwmethod geometrized_to_CGS(magnitude, dim, configurations::Configurations;) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=configurations.unit_mass_in_solar_masses)
@kwmethod CGS_to_geometrized(magnitude, dim, configurations::Configurations;) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=configurations.unit_mass_in_solar_masses)

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