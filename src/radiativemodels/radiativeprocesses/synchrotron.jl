function synchrotron_M(x_M, α, β, γ)
    return (4.0505 * α / x_M^(1/6)) * (1 + (0.40 * β / x_M^(1/4)) + (0.5316 * γ / x_M^(1/2))) * exp(-1.8899 * x_M^(1/3))
end

function synchrotron_x_M(ν, ν0, θe)
    return 2ν/(3ν0*θe^2)
end

function critical_frequency(ν0, ne, θe, B, K2, R)
    f_νC(xM) = exp(1.8899 * xM^(1 / 3)) - 2.49e-10 * 4.0 * π * ne * R * (1.0 / xM^(1 / 6) + 0.4 / xM^(17 / 12) + 0.5316 / xM^(5 / 3)) / (B * θe^3 * K2)
    xM = find_zero(f_νc, (1.0e4, 1.0e12), Bisection())
    νc = 1.5 * ν0 * θe^2 * xM
    return νc
end
#TODO evaluate black body below critical frequency
function synchrotron_emissivity(ε, ne, Te, B, α, β, γ)
    h = PhysicalConstants.h
    k_B = PhysicalConstants.k_B
    ce = PhysicalConstants.ce
    c = PhysicalConstants.c
    me = PhysicalConstants.me
    ν = ε/h
    ν0 = (ce*B)/(2π*me*c)
    θe = k_B*Te/(me*c^2)
    x_M = synchrotron_x_M(ν, ν0, θe)
    if Te > 3.2e10
        K2 = besselk(2,1/θe)
        # νc = critical_frequency(ν0, ne, θe, B, K2, R)
        jε = (ce^2/(h*c*sqrt(3)))*(ne*ν/K2)*synchrotron_M(x_M, 1.0, 1.0, 1.0)
    elseif Te > 5e8
        K2 = besselk(2,1/θe)
        # νc = critical_frequency(ν0, ne, θe, B, K2, R)
        jε = (ce^2/(h*c*sqrt(3)))*(ne*ν/K2)*synchrotron_M(x_M, α, β, γ)
    else
        jε = 0
    end
    return jε
end