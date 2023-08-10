module PhysicalConstants

#Physical and astrophysical constants
# Values taken from Particle Data Group "Physical constants" and "Astrophysical constants"
# https://pdg.lbl.gov/  (last updated  2021)

const c = 2.99792458e10         #Speed of light in CGS                      (exact)
const c2 = 8.98755178737e20     #Squared speed of light in CGS              (exact)    
const h = 6.62607015e-27        #Planck constant in CGS  6.626 070 15       (exact)
const h_eV = 4.13566769692e-15  #Planck's constant in eV*s                  (exact*)
const h_bar = 1.054571817e-27   #Reduced Planck constant in CGS             (exact*) 
const k_B = 1.380649e-16        #Boltzmann constant in CGS                  (exact)    
const σ = 5.670374419e-5        #Stefan-Boltzmann constant in CGS           (exact*)
const eV = 1.602176634e-12      #1eV in erg                                 (exact)
const me = 9.109383701528e-28   #Electron mass in CGS 
const ce = 4.80320471257e-10    #Electron charge in esu-CGS and Gaussian-CGS (statC)   (exact*) 
#const ce = 1.602176634e-20     #Electron charge in emu-CGS (a relic, just in case)    (exact)
const G = 6.6743015e-8          #Newton gravitational constant in CGS
const M_sun = 1.988414e33       #Solar mass in CGS
const pc = 3.08567758149e18     #1pc in cm                                             (exact)
const mu = 1.6605390666050e-24  #Atomic mass unit in CGS
end