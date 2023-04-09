export geometrized_to_CGS
export CGS_to_geometrized

geometrized_to_CGS(magnitude; dim, M1::Number) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=M1) 
CGS_to_geometrized(magnitude; dim, M1::Number) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=M1)

geometrized_to_CGS(magnitude; dim, configurations::Configurations) = magnitude*geometrized_unit_magnitude_in_CGS(dim; M1=configurations.M1) 
CGS_to_geometrized(magnitude; dim, configurations::Configurations) = magnitude/geometrized_unit_magnitude_in_CGS(dim; M1=configurations.M1)

function geometrized_unit_magnitude_in_CGS(dim; M1)

    #Returns the value in CGS of the geometrized system unit of dimension dim
    #M1 is the chosen unit mass in solar masses 

    c = PhysicalConstants.c
    G = PhysicalConstants.G
    M_sun = PhysicalConstants.M_sun

    M1_cgs = M1*M_sun        #Geometrized unit mass in CGS 
    L1_cgs = G*M1_cgs/c^2    #Geometrized unit length in CGS
    T1_cgs = L1_CGS/c        #Geometrized unit time in CGS

    return L1_cgs^dim.L*M1_cgs^dim.M*T1_cgs^dim.T

end