function synchrotron_M(x_M, α, β, γ)
    return (4.0505 * α / x_M^(1 / 6)) *
           (1 + (0.40 * β / x_M^(1 / 4)) + (0.5316 * γ / x_M^(1 / 2))) *
           exp(-1.8899 * x_M^(1 / 3))
end

function synchrotron_x_M(ν, ν0, θe)
    return 2ν / (3ν0 * θe^2)
end

function critical_frequency(ν0, ne, θe, B, K2, R)
    all(x -> isfinite(x) && x>0,(ν0,ne,θe,B,K2,R)) ||
        throw(ArgumentError("Critical-frequency inputs must be finite and positive."))
    logC = log(2.49e-10*4π)+log(ne)+log(R)-log(B)-3log(θe)-log(K2)
    # Solve in log xM: the old exponential overflowed before root finding.
    # Dividing the optically thin emissivity by the Rayleigh-Jeans intensity
    # gives xM^(-7/6), not xM^(-1/6), in the first term.
    f(logx) = 1.8899*exp(logx/3)-logC+7logx/6-
        log1p(0.4exp(-logx/4)+0.5316exp(-logx/2))
    xM = exp(find_zero(f,(-60.0,60.0),Bisection()))
    νc = 1.5 * ν0 * θe^2 * xM
    return νc
end
# Straub/Mahadevan fit: below 5e8 K the legacy approximation neglects
# cyclo-synchrotron emission. It is a documented cutoff, not an extrapolation.
function synchrotron_coefficients(sy::Union{Synchrotron,SynchrotronAndBremsstrahlung},Te)
    Te < 5e8 && return (0.0,0.0,0.0)
    Te >= 3.2e10 && return (1.0,1.0,1.0)
    if Te > last(sy.T)
        weight = (Te-last(sy.T))/(3.2e10-last(sy.T))
        return ((1-weight)*sy.α(last(sy.T))+weight,
            (1-weight)*sy.β(last(sy.T))+weight,
            (1-weight)*sy.γ(last(sy.T))+weight)
    end
    return sy.α(Te),sy.β(Te),sy.γ(Te)
end

function synchrotron_emissivity(ε, ne, Te, B, α, β, γ)
    ε > 0 || throw(DomainError(ε,"Photon energy must be positive."))
    ne >= 0 && Te >= 0 && B >= 0 || throw(ArgumentError("Density, temperature and magnetic field must be nonnegative."))
    (Te < 5e8 || ne == 0 || B == 0) && return 0.0
    h = PhysicalConstants.h
    k_B = PhysicalConstants.k_B
    ce = PhysicalConstants.ce
    c = PhysicalConstants.c
    me = PhysicalConstants.me
    ν = ε / h
    ν0 = (ce * B) / (2π * me * c)
    θe = k_B * Te / (me * c^2)
    x_M = synchrotron_x_M(ν, ν0, θe)
    Te >= 3.2e10 && ((α,β,γ) = (1.0,1.0,1.0))
    K2 = besselk(2, 1 / θe)
    return (ce^2 / (h * c * sqrt(3))) * (ne * ν / K2) * synchrotron_M(x_M, α, β, γ)
end

function synchrotron_emissivity(ε,ne,Te,B,process::Union{Synchrotron,SynchrotronAndBremsstrahlung})
    coefficients = synchrotron_coefficients(process,Te)
    return synchrotron_emissivity(ε,ne,Te,B,coefficients...)
end
