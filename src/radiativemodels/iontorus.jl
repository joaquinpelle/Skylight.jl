"""
Ion torus model from https://www.aanda.org/articles/aa/abs/2012/07/aa19209-12/aa19209-12.html
"""

@with_kw mutable struct IonTorus{T} <: AbstractRadiativeModel
    λ::Float64 = 0.7 #Specific angular momentum dimensionless parameter
    εc::Float64 = 1e-17 #Central energy density in CGS
    n::Float64 = 3/2 #Polytropic index
    Tec::Float64 = 1e10 #Central electorn temperature in Kelvin
    ξ::Float64 = 0.01 #Electron to proton temperature ratio at the center 
    H_abundance::Float64 = 0.75
    He_abundance::Float64 = 0.25
    β::Float64 = 0.45 #Equipartition factor
    rotation_sense::T = ProgradeRotation()
    μᵢ::Float64 = 4/(4H_abundance+He_abundance)
    μₑ::Float64 = 2/(1+He_abundance)
    𝓜₀::Float64 = μᵢ/(μᵢ+μₑ)
    𝓜₁::Float64 = μᵢ*ξ/(μᵢ*ξ+μₑ)
    K::Float64 = PhysicalConstants.k_B*Tec/((1-β)*PhysicalConstants.mu*εc^(1/n)*μₑ*𝓜₁) 
    Hc::Float64 = (n+1)*log(1+K*εc^(1/n))
    l0::Float64 = 0.0
    rcusp::Float64 = 0.0
    rcenter::Float64 = 0.0
    potential_at_surface::Float64 = 0.0
    potential_at_center::Float64 = 0.0
    _l0_is_set::Bool = false
    _radii_are_set::Bool = false
    @assert 0 ≤ λ ≤ 1 "λ must be between 0 and 1"
    @assert isa(rotation_sense, AbstractRotationSense) "Rotation sense must be either ProgradeRotation() or RetrogradeRotation()"
end

function torus_specific_angular_momentum!(model::IonTorus, spacetime)
    λ = model.λ
    lms = innermost_stable_circular_orbit_specific_angular_momentum(spacetime, model.rotation_sense)
    lmb = marginally_bound_specific_angular_momentum(spacetime, model.rotation_sense)
    model.l0 = lms + λ*(lmb - lms)
    model._l0_is_set = true
    return nothing
end

function cusp_and_center_radius!(model::IonTorus, spacetime)
    if !model._l0_is_set
        println("Specific angular momentum not set")
        return nothing
    end
    M = spacetime.M
    l0 = model.l0
    rmax = 100M
    aux(r) = circular_geodesic_specific_angular_momentum([0.0,r,π/2,0.0], spacetime, model.rotation_sense)-l0
    model.rcusp = find_zero(aux, (rmb, rms))
    model.rmax = find_zero(aux, (rms, rmax))
    model._radii_are_set = true
    return nothing
end

function constant_angular_momentum_angular_speed(position, spacetime, model::IonTorus, g)
    l0 = model.l0
    metric!(g, position, spacetime)
    gtt = g[1,1]
    gtφ = g[1,4]
    gφφ = g[4,4]
    return -(gtφ+l0*gtt)/(gφφ+l0*gtφ)
end 

function torus_potential(position, spacetime, model::IonTorus, g)
    l0 = model.l0
    metric!(g, position, spacetime)
    gtt = g[1,1]
    gtφ = g[1,4]
    gφφ = g[4,4]
    Ω = -(gtφ+l0*gtt)/(gφφ+l0*gtφ)
    return 0.5*log(-(gtt+2Ω*gtφ+Ω^2*gφφ)/(gtt+Ω*gtφ)^2)
end

torus_potential(position, spacetime, model::IonTorus) = torus_potential(position, spacetime, model, zeros(4,4))

function torus_potential_at_surface(spacetime, model::IonTorus)
    position = equatorial_position(model.rcusp, spacetime) 
    return torus_potential(position, spacetime, model) 
end

function torus_potential_at_center(spacetime, model::IonTorus)
    return model.Hc+torus_potential_at_surface(spacetime, model) 
end

function torus_potentials_at_center_and_surface!(model::IonTorus, spacetime)
    model.potential_at_surface = torus_potential_at_surface(spacetime, model)
    model.potential_at_center = torus_potential_at_center(spacetime, model)
    return nothing
end

function torus_normalized_potential(position, spacetime, model::IonTorus, g)
    Ws = model.potential_at_surface
    Wc = model.potential_at_center
    W = torus_potential(position, spacetime, model, g)
    return (W-Ws)/(Wc-Ws)
end

#TODO evaluate omega>0 conditional
function energy_density(position, spacetime, model::IonTorus, g)
    n = model.n
    K = model.K
    εc = model.εc
    ω = torus_normalized_potential(position, spacetime, model, g)
    ε = K^(-n)*((K*εc^(1/n)+1)^ω - 1)^n
    return ε
end

function pressure(position, spacetime, model::IonTorus, g)
    n = model.n
    K = model.K
    ε = energy_density(position, spacetime, model, g)
    return K*ε^(1+1/n)
end

function electron_temperature(position, spacetime, model::IonTorus, g)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    n = model.n
    K = model.K
    εc = model.εc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    ω = torus_normalized_potential(position, spacetime, model, g)
    ε = K^(-n)*((K*εc^(1/n)+1)^ω - 1)^n
    P = K*ε^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ε)
    return ((1-ω)*𝓜₀+ω*𝓜₁)*μₑ*factor
end

function ion_temperature(position, spacetime, model::IonTorus, g)
    K = model.K
    n = model.polytropic_index
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    n = model.n
    K = model.K
    εc = model.εc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    ω = torus_normalized_potential(position, spacetime, model, g)
    ε = K^(-n)*((K*εc^(1/n)+1)^ω - 1)^n
    P = K*ε^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ε)
    return ((μₑ/μᵢ)𝓜₀+ω*(𝓜₀-𝓜₁))*μᵢ*factor
end

function dimensionless_electron_temperature(position, spacetime, model::IonTorus, g)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    c = PhysicalConstants.c
    me = PhysicalConstants.me
    n = model.n
    K = model.K
    εc = model.εc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    ω = torus_normalized_potential(position, spacetime, model, g)
    ε = K^(-n)*((K*εc^(1/n)+1)^ω - 1)^n
    P = K*ε^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ε)
    Te = ((1-ω)*𝓜₀+ω*𝓜₁)*μₑ*factor
    return k_B*Te/(me*c^2)
end

function dimensionless_ion_temperature(position, spacetime, model::IonTorus, g)
    K = model.K
    n = model.polytropic_index
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    c = PhysicalConstants.c
    n = model.n
    K = model.K
    εc = model.εc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    mi = mu*μᵢ
    ω = torus_normalized_potential(position, spacetime, model, g)
    ε = K^(-n)*((K*εc^(1/n)+1)^ω - 1)^n
    P = K*ε^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ε)
    Ti = ((μₑ/μᵢ)𝓜₀+ω*(𝓜₀-𝓜₁))*μᵢ*factor
    return k_B*Ti/(mi*c^2)
end

torus_normalized_potential(position, spacetime, model::IonTorus) = torus_normalized_potential(position, spacetime, model, zeros(4,4))
energy_density(position, spacetime, model::IonTorus) = energy_density(position, spacetime, model, zeros(4,4))
pressure(position, spacetime, model::IonTorus) = pressure(position, spacetime, model, zeros(4,4))
electron_temperature(position, spacetime, model::IonTorus) = electron_temperature(position, spacetime, model, zeros(4,4))
ion_temperature(position, spacetime, model::IonTorus) = ion_temperature(position, spacetime, model, zeros(4,4))
dimensionless_electron_temperature(position, spacetime, model::IonTorus) = dimensionless_electron_temperature(position, spacetime, model, zeros(4,4))
dimensionless_ion_temperature(position, spacetime, model::IonTorus) = dimensionless_ion_temperature(position, spacetime, model, zeros(4,4))