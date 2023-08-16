"""
Ion torus model from https://www.aanda.org/articles/aa/abs/2012/07/aa19209-12/aa19209-12.html
"""

@with_kw mutable struct IonTorus{R,P} <: AbstractRadiativeModel
    λ::Float64 = 0.3 #Specific angular momentum dimensionless parameter
    ϵc::Float64 = 1e-17 #Central energy density in CGS
    n::Float64 = 3/2 #Polytropic index
    Tec::Float64 = 2e9 #Central electorn temperature in Kelvin
    ξ::Float64 = 0.1 #Electron to proton temperature ratio at the center 
    β::Float64 = 0.1 #Equipartition factor
    H_abundance::Float64 = 0.75
    He_abundance::Float64 = 0.25
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

function IonTorus(spacetime::AbstractSpacetime; kwargs...)
    model = IonTorus(; kwargs...)
    torus_specific_angular_momentum!(model, spacetime)
    cusp_and_center_radius!(model, spacetime)
    torus_potentials_at_center_and_surface!(model, spacetime)
    return model
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

function constant_angular_momentum_angular_speed(g::AbstractMatrix, model::IonTorus)
    l0 = model.l0
    gtt = g[1,1]
    gtφ = g[1,4]
    gφφ = g[4,4]
    return -(gtφ+l0*gtt)/(gφφ+l0*gtφ)
end 

function torus_potential(g::AbstractMatrix, model::IonTorus)
    gtt = g[1,1]
    gtφ = g[1,4]
    gφφ = g[4,4]
    Ω = constant_angular_momentum_angular_speed(g, model)
    p2 = gtt+2Ω*gtφ+Ω^2*gφφ
    return ifelse(p2<0.0, 0.5*log(abs(p2)/(gtt+Ω*gtφ)^2), -Inf)
end

function torus_potential_at_surface(spacetime, model::IonTorus)
    position = equatorial_position(model.rcusp, coordinates_topology(spacetime)) 
    g = metric(position, spacetime)
    return torus_potential(g, model) 
end

function torus_potential_at_center(spacetime, model::IonTorus)
    # return model.Hc+torus_potential_at_surface(spacetime, model) 
    position = equatorial_position(model.rcenter, coordinates_topology(spacetime)) 
    g = metric(position, spacetime)
    return torus_potential(g, model) 
end

function torus_potentials_at_center_and_surface!(model::IonTorus, spacetime)
    model.potential_at_surface = torus_potential_at_surface(spacetime, model)
    model.potential_at_center = torus_potential_at_center(spacetime, model)
    return nothing
end

function torus_normalized_potential(W::Real, model::IonTorus)
    Ws = model.potential_at_surface
    Wc = model.potential_at_center
    ω = (W-Ws)/(Wc-Ws)
    return ω
end

function energy_density(ω::Real, model::IonTorus)
    n = model.n
    K = model.K
    ϵc = model.ϵc
    ϵ = ifelse(ω>0.0, K^(-n)*((K*ϵc^(1/n)+1)^ω - 1)^n, 0.0)
    return ϵ
end

function pressure(ϵ::Real, model::IonTorus)
    n = model.n
    K = model.K
    return K*ϵ^(1+1/n)
end

function electron_temperature(ω::Real, ϵ::Real, model::IonTorus)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    K = model.K
    n = model.n
    factor = (1-β)*mu*K*ϵ^(1/n)/k_B
    return ((1-ω)*𝓜₀+ω*𝓜₁)*μₑ*factor
end

function ion_temperature(ω::Real, ϵ::Real, model::IonTorus)
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    K = model.K
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    factor = (1-β)*mu*K*ϵ^(1/n)/k_B
    return ((μₑ/μᵢ)𝓜₀+ω*(𝓜₀-𝓜₁))*μᵢ*factor
end

function number_densities(ϵ::Real, model::IonTorus)
    mu = PhysicalConstants.mu
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    ne = ϵ/(mu*μₑ)
    ni = ϵ/(mu*μᵢ)
    return ne, ni
end

function electron_number_density(ϵ::Real, model::IonTorus)
    mu = PhysicalConstants.mu
    μₑ = model.μₑ
    ne = ϵ/(mu*μₑ)
    return ne
end

function magnetic_field(P::Real, model::IonTorus)
    β = model.β
    return sqrt(24π*β*P)
end

function number_densities_and_electron_temperature(ω::Real, model::IonTorus)
    ϵ = energy_density(ω, model)
    Te = electron_temperature(ω, ϵ, model)
    ne, ni = number_densities(ϵ, model)
    return ne, ni, Te
end

function electron_number_density_temperature_and_magnetic_field(ω, model::IonTorus)
    ϵ = energy_density(ω, model)
    Te = electron_temperature(ω, ϵ, model)
    ne = electron_number_density(ϵ, model)
    B = magnetic_field(P, model)
    return ne, Te, B
end

#TODO evaluate doing inside source in equations.jl
function rest_frame_four_velocity!(vector, position, metric, spacetime, model::IonTorus, coords_top)
    Ω = constant_angular_momentum_angular_speed(metric, model)
    p2 = metric[1,1]+2Ω*metric[1,4]+Ω^2*metric[4,4]
    ifelse(p2 < 0.0,
        circular_motion_four_velocity_allowing_spacelike!(vector, position, Ω, metric, coords_top),
        static_four_velocity!(vector, metric))
end

rest_frame_absorptivity!(αε, position, ε, g, spacetime, model::IonTorus, coords_top) = nothing
rest_frame_emissivity!(jε, position, ε, g, spacetime, model::IonTorus, coords_top) = rest_frame_emissivity!(model.radiative_process, jε, position, ε, g, spacetime, model, coords_top)

function rest_frame_emissivity!(::Bremsstrahlung, jε, position, ε, g, spacetime, model::IonTorus, coords_top)
    ω = torus_normalized_potential(g, model)
    if ω>0
        ne, ni, Te = number_densities_and_electron_temperature(ω, model)
        for (i,εk) in enumerate(ε)
            jε[i] = bremsstrahlung_emissivity(εk, ne, ni, Te)
        end
    end
    return nothing
end

#TODO revise that I probably should not calculate the variables outside the torus
# function rest_frame_emissivity!(sy::Synchrotron, jε, position, ε, g, spacetime, model::IonTorus, coords_top)
#     ω = torus_normalized_potential(position, spacetime, model, g)
#     ne, Te, B = electron_number_density_temperature_and_magentic_field(model, ω)
#     α = sy.α(Te)
#     β = sy.β(Te)
#     γ = sy.γ(Te)
#     for (i,εk) in enumerate(ε)
#         jε[i] = ifelse(ω>0, 
#                        synchrotron_emissivity(εk, ne, Te, B, α, β, γ),
#                        0.0)
#     end
#     return nothing
# end

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
function torus_normalized_potential(g::AbstractMatrix, model::IonTorus)
    W = torus_potential(g, model)
    return torus_normalized_potential(W, model)
end

function energy_density(g::AbstractMatrix, model::IonTorus)
    ω = torus_normalized_potential(g, model)
    return energy_density(ω, model)
end

function pressure(g::AbstractMatrix, model::IonTorus)
    ϵ = energy_density(g, model)
    return pressure(ϵ, model)
end

function electron_temperature(g::AbstractMatrix, model::IonTorus)
    ω = torus_normalized_potential(g, model)
    ϵ = energy_density(ω, model) 
    return electron_temperature(ω, ϵ, model)
end

function ion_temperature(g::AbstractMatrix, model::IonTorus)
    ω = torus_normalized_potential(g, model)
    ϵ = energy_density(ω, model) 
    P = pressure(ϵ, model)
    return ion_temperature(ω, ϵ, model)
end

function number_densities(g::AbstractMatrix, model::IonTorus)
    ϵ = energy_density(g, model)
    return number_densities(ϵ, model)
end

torus_potential(position, spacetime, model::IonTorus) = torus_potential(metric(position, spacetime), model)
torus_normalized_potential(position, spacetime, model::IonTorus) = torus_normalized_potential(metric(position, spacetime), model)
energy_density(position, spacetime, model::IonTorus) = energy_density(metric(position, spacetime), model)
pressure(position, spacetime, model::IonTorus) = pressure(metric(position, spacetime), model)
electron_temperature(position, spacetime, model::IonTorus) = electron_temperature(metric(position, spacetime), model) 
ion_temperature(position, spacetime, model::IonTorus) = ion_temperature(metric(position, spacetime), model) 
number_densities(position, spacetime, model::IonTorus) = number_densities(metric(position, spacetime), model)
magnetic_field(position, spacetime, model::IonTorus) = magnetic_field(pressure(position, spacetime, model), model)
number_densities_and_electron_temperature(position, spacetime, model::IonTorus) = number_densities_and_electron_temperature(torus_normalized_potential(position, spacetime, model), model)
electron_number_density_temperature_and_magnetic_field(position, spacetime, model::IonTorus) = electron_number_density_temperature_and_magnetic_field(torus_normalized_potential(position, spacetime, model), model)