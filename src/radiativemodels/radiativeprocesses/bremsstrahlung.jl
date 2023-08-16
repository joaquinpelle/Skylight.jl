"""CGS units. Functions taken from Straub et al (2012), but jε instead of jν"""

function bremsstrahlung_Fei(θe)
    return ifelse(θe < 1, 
                  4(2θe/π^3)^0.5*(1+1.781θe^1.34),
                  9*θe/(2π)*(log(1.123θe+0.48)+1.5))
end

function bremsstrahlung_Fee(θe)
    η = 0.5616
    return ifelse(θe < 1, 
                  20/(9*sqrt(π))*(44-3π^2*θe^1.5*(1+1.1*θe+θe^2-1.25*θe^2.5)),
                  24*θe*(log(2*η*θe)+1.28))
end

function bremsstrahlung_fee(ne, θe)
    re = PhysicalConstants.re
    c = PhysicalConstants.c
    αf = PhysicalConstants.alpha_f
    me = PhysicalConstants.me
    return ne^2*re^2*αf*me*c^3*bremsstrahlung_Fee(θe)
end

function bremsstrahlung_fei(ne, ni, θe)
    σT = PhysicalConstants.sigma_T
    c = PhysicalConstants.c
    αf = PhysicalConstants.alpha_f
    me = PhysicalConstants.me
    return ne*ni*σT*αf*me*c^3*bremsstrahlung_Fei(θe)
end

function bremsstrahlung_emissivity(ε, ne, ni, Te)
    γE = Base.MathConstants.eulergamma 
    k_B = PhysicalConstants.k_B
    c = PhysicalConstants.c
    me = PhysicalConstants.me
    θe = k_B*Te/(me*c^2)
    x = k_B*Te/ε
    try
        Gmean = ifelse(x < 1, 
                    3*x/π,
                    sqrt(3)/π*log(4/γE*x))
        fbr = bremsstrahlung_fee(ne, θe) + bremsstrahlung_fei(ne, ni, θe)
        jε = fbr/(4π*k_B*Te)*exp(-1/x)*Gmean
    return jε 
    catch e
        println("Te = $Te", " x = $x")
        rethrow(e)
    end
end