"""
Ion torus model from https://www.aanda.org/articles/aa/abs/2012/07/aa19209-12/aa19209-12.html
"""

@with_kw mutable struct IonTorus{R,P} <: AbstractRadiativeModel
    λ::Float64 = 0.7 #Specific angular momentum dimensionless parameter
    ϵc::Float64 = 1e-17 #Central energy density in CGS
    n::Float64 = 3/2 #Polytropic index
    Tec::Float64 = 1e10 #Central electorn temperature in Kelvin
    ξ::Float64 = 0.01 #Electron to proton temperature ratio at the center 
    H_abundance::Float64 = 0.75
    He_abundance::Float64 = 0.25
    β::Float64 = 0.45 #Equipartition factor
    rotation_sense::R = ProgradeRotation()
    radiative_process::P = Bremsstrahlung()
    μᵢ::Float64 = 4/(4H_abundance+He_abundance)
    μₑ::Float64 = 2/(1+He_abundance)
    𝓜₀::Float64 = μᵢ/(μᵢ+μₑ)
    𝓜₁::Float64 = μᵢ*ξ/(μᵢ*ξ+μₑ)
    K::Float64 = PhysicalConstants.k_B*Tec/((1-β)*PhysicalConstants.mu*ϵc^(1/n)*μₑ*𝓜₁) 
    Hc::Float64 = (n+1)*log(1+K*ϵc^(1/n))
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
    lms = innermost_stable_specific_angular_momentum(spacetime, model.rotation_sense)
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
    rmb = mbco_radius(spacetime, model.rotation_sense)
    rms = isco_radius(spacetime, model.rotation_sense)
    rmax = 100M
    aux(r) = circular_geodesic_specific_angular_momentum([0.0,r,π/2,0.0], spacetime, model.rotation_sense)-l0
    model.rcusp = find_zero(aux, (rmb, rms))
    model.rcenter = find_zero(aux, (rms, rmax))
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
    position = equatorial_position(model.rcusp, coordinates_topology(spacetime)) 
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

#TODO evaluate omega>0 conditional with ifelse()
function energy_density(position, spacetime, model::IonTorus, g)
    n = model.n
    K = model.K
    ϵc = model.ϵc
    ω = torus_normalized_potential(position, spacetime, model, g)
    ϵ = K^(-n)*((K*ϵc^(1/n)+1)^ω - 1)^n
    return ϵ
end

function pressure(position, spacetime, model::IonTorus, g)
    n = model.n
    K = model.K
    ϵ = energy_density(position, spacetime, model, g)
    return K*ϵ^(1+1/n)
end

function electron_temperature(position, spacetime, model::IonTorus, g)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    n = model.n
    K = model.K
    ϵc = model.ϵc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    ω = torus_normalized_potential(position, spacetime, model, g)
    ϵ = K^(-n)*((K*ϵc^(1/n)+1)^ω - 1)^n
    P = K*ϵ^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ϵ)
    return ((1-ω)*𝓜₀+ω*𝓜₁)*μₑ*factor
end

function ion_temperature(position, spacetime, model::IonTorus, g)
    K = model.K
    n = model.polytropic_index
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    n = model.n
    K = model.K
    ϵc = model.ϵc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    ω = torus_normalized_potential(position, spacetime, model, g)
    ϵ = K^(-n)*((K*ϵc^(1/n)+1)^ω - 1)^n
    P = K*ϵ^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ϵ)
    return ((μₑ/μᵢ)𝓜₀+ω*(𝓜₀-𝓜₁))*μᵢ*factor
end

function number_densities(position, spacetime, model::IonTorus, g)
    mu = PhysicalConstants.mu
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    ϵ = energy_density(position, spacetime, model, g)
    ne = ϵ/(mu*μₑ)
    ni = ϵ/(mu*μᵢ)
    return ne, ni
end

function number_densities_and_electron_temperature(model::IonTorus, ω)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    n = model.n
    K = model.K
    ϵc = model.ϵc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    ϵ = K^(-n)*((K*ϵc^(1/n)+1)^ω - 1)^n
    P = K*ϵ^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ϵ)
    ne = ϵ/(mu*μₑ)
    ni = ϵ/(mu*μᵢ)
    Te = ((1-ω)*𝓜₀+ω*𝓜₁)*μₑ*factor
    return ne, ni, Te
end

function electron_number_density_temperature_and_magnetic_field(model::IonTorus, ω)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    μₑ = model.μₑ
    n = model.n
    K = model.K
    ϵc = model.ϵc
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    ϵ = K^(-n)*((K*ϵc^(1/n)+1)^ω - 1)^n
    P = K*ϵ^(1+1/n)
    factor = (1-β)*mu*P/(k_B*ϵ)
    ne = ϵ/(mu*μₑ)
    Te = ((1-ω)*𝓜₀+ω*𝓜₁)*μₑ*factor
    B = sqrt(24π*β*P)
    return ne, Te, B
end

#TODO rename rest_frame functions
#TODO beware superluminal four v
function rest_frame_four_velocity!(vector, position, metric, spacetime, model::IonTorus, coords_top)
    angular_speed = constant_angular_momentum_angular_speed(position, spacetime, model, metric)
    circular_motion_four_velocity!(vector, position, angular_speed, metric, coords_top)
end
rest_frame_absorptivity!(αε, position, ε, g, spacetime, model::IonTorus, coords_top) = nothing
rest_frame_emissivity!(jε, position, ε, g, spacetime, model::IonTorus, coords_top) = rest_frame_emissivity!(model.radiative_process, jε, position, ε, g, spacetime, model, coords_top)
#TODO benchmark ifelse against regular if
function rest_frame_emissivity!(::Bremsstrahlung, jε, position, ε, g, spacetime, model::IonTorus, coords_top)
    ω = torus_normalized_potential(position, spacetime, model, g)
    ne, ni, Te = number_densities_and_electron_temperature(model, ω)
    for (i,εk) in enumerate(ε)
        jε[i] = ifelse(ω>0, 
                       bremsstrahlung_emissivity(εk, ne, ni, Te),
                       0.0)
    end
    return nothing
end

#TODO revise that I probably should not calculate the variables outside the torus
function rest_frame_emissivity!(sy::Synchrotron, jε, position, ε, g, spacetime, model::IonTorus, coords_top)
    ω = torus_normalized_potential(position, spacetime, model, g)
    ne, Te, B = electron_number_density_temperature_and_magentic_field(model, ω)
    α = sy.α(Te)
    β = sy.β(Te)
    γ = sy.γ(Te)
    for (i,εk) in enumerate(ε)
        jε[i] = ifelse(ω>0, 
                       synchrotron_emissivity(εk, ne, Te, B, α, β, γ),
                       0.0)
    end
    return nothing
end

#TODO: implement total
# function rest_frame_emissivity!(::Total, jε, position, ε, g, spacetime, model::IonTorus, coords_top)
#     ω = torus_normalized_potential(position, spacetime, model, g)
#     ne, ni, Te = electron_number_density_temperature_and_magentic_field(model, ω)
#     for (i,εk) in enumerate(ε)
#         jε[i] = ifelse(ω>0, 
#                        bremsstrahlung_emissivity(ne, ni, Te, εk),
#                        0.0)
#     end
#     return nothing
# end
torus_normalized_potential(position, spacetime, model::IonTorus) = torus_normalized_potential(position, spacetime, model, zeros(4,4))
energy_density(position, spacetime, model::IonTorus) = energy_density(position, spacetime, model, zeros(4,4))
pressure(position, spacetime, model::IonTorus) = pressure(position, spacetime, model, zeros(4,4))
electron_temperature(position, spacetime, model::IonTorus) = electron_temperature(position, spacetime, model, zeros(4,4))
ion_temperature(position, spacetime, model::IonTorus) = ion_temperature(position, spacetime, model, zeros(4,4))
number_densities(position, spacetime, model::IonTorus) = number_densities(position, spacetime, model, zeros(4,4))