function synchrotron_M(x_M, α, β, γ)
    return (4.0505 * α / x_M^(1/6)) * (1 + (0.40 * β / x_M^(1/4)) + (0.5316 * γ / x_M^(1/2))) * exp(-1.8899 * x_M^(1/3))
end

function synchrotron_x_M(ν, ν0, θe)
    return 2ν / (3ν0 * θe^2)
end

function synchrotron_emissivity(ε, ne, Te, B, α, β, γ)
    h = PhysicalConstants.h
    k_B = PhysicalConstants.k_B
    e = PhysicalConstants.e
    c = PhysicalConstants.c
    me = PhysicalConstants.me
    ν = ε/h
    ν0 = (e*B)/(2π*me*c)
    θe = k_B*Te/(me*c^2)
    x_M = synchrotron_x_M(ν, ν0, θe)
    return (e^2/(h*c*√3))*(ne*ν/besselk(2,1/θe))*M(x_M, α, β, γ)
end