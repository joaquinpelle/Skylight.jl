thermal_emission_bolometric_intensity(T) = stefan_boltzmann_function(T)/π
thermal_emission_window_intensity(E1, E2, T) = planck_integral(E1, T) - planck_integral(E2, T)
thermal_emission_specific_intensity(E, T) = planck_function(E, T)

function stefan_boltzmann_function(T)
    σ = PhysicalConstants.σ
    return σ*T^4
end

"""
    planck_function(E, T)

Planck's law of black body radiation in terms of energy.

# Arguments
- `E::Float64`: The energy in erg.
- `T::Float64`: The temperature in Kelvin.

# Returns
- `::Float64`: The planck function in CGS.
"""
function planck_function(E,T)
    c = PhysicalConstants.c      
    h = PhysicalConstants.h      
    k_B = PhysicalConstants.k_B  
    return 2*E^3/(h^3*c^2*expm1(E/(k_B*T)))
end

"""
    planck_integral(E, T)

Estimation of the Planck integral from E to infinity for temperature T. 
Based on "Planck functions and integrals; methods of computation" (T.E. Michels)
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680008986.pdf

# Arguments
- `E::Float64`: The energy in erg.
- `T::Float64`: The temperature in Kelvin.

# Returns
- `::Float64`: The planck integral in CGS.

# Remarks
An evaluation of the error is too expensive for millions of calls, so we fix N=5.
This was checked against scipy integration, yielding at least 3 siginficant digits for all x.
When x>1 the amount of siginificant digits is 8 or more
"""
function planck_integral(E,T)
    c = PhysicalConstants.c     
    h = PhysicalConstants.h     
    k_B = PhysicalConstants.k_B 
    
    factor = 2*(k_B*T)^4/(h^3*c^2)  
    
    x = E/(k_B*T) 
    x2 = x*x
    x3 = x*x2
 
    return  factor*sum(exp(-n*x)*(x3 + (3*x2 + 6*(x+1/n)/n)/n)/n for n in 1:4)
end

"""
    rayleigh_jeans_function(E, T)

Rayleigh-Jeans law of black body radiation in terms of energy.

# Arguments
- `E::Float64`: The energy in erg.
- `T::Float64`: The temperature in Kelvin.

# Returns
- `::Float64`: The Rayleigh-Jeans function in CGS.
"""
function rayleigh_jeans_function(E, T)
    c = PhysicalConstants.c
    h = PhysicalConstants.h
    k_B = PhysicalConstants.k_B
    return (8 * π / (c^3 * h^3)) * E^2 * k_B * T
end