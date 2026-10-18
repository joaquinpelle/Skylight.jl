"""
    IonTorus <: AbstractRadiativeModel

Geometrically thick, optically thin ion torus model with synchrotron and bremsstahlung emission ([Straub et al. 2012](https://www.aanda.org/articles/aa/abs/2012/07/aa19209-12/aa19209-12.html)).

# Fields
- `λ::Float64`: Specific angular momentum dimensionless parameter
- `ϵc::Float64`: Central mass density in g cm⁻³ (legacy field name)
- `n::Float64`: Polytropic index
- `Tec::Float64`: Central electron temperature in Kelvin
- `ξ::Float64`: Electron to proton temperature ratio at the center
- `β::Float64`: Equipartition factor
- `H_abundance::Float64`: Hydrogen abundance
- `He_abundance::Float64`: Helium abundance
- `rotation_sense::R`: Sense of rotation of the torus, which can be either `ProgradeRotation()` or `RetrogradeRotation()`
- `radiative_process::P`: Radiative process, which can be either `Bremsstrahlung()`, `Synchrotron()` or `SynchrotronAndBremsstrahlung()`

# Constructor
```julia
spacetime = KerrSpacetimeBoyerLindquistCoordinates(M = 1.0, a = 0.7)
IonTorus = IonTorus(spacetime; 
    λ = 0.3, 
    ϵc = 1e-17, 
    n = 3 / 2, 
    ξ = 0.1, 
    β = 0.1, 
    H_abundance = 0.75, 
    He_abundance = 0.25, 
    rotation_sense = ProgradeRotation(), 
    radiative_process = Bremsstrahlung())
```

# Thermodynamics
`IonTorus(spacetime)` derives `Tec` from hydrostatic balance, using e ≈ ρc² and
a total-pressure polytrope. An explicitly supplied `Tec` must agree with that
balance. Use `thermodynamics=:prescribed, Tec=...` only for a phenomenological
temperature profile that does not in general satisfy hydrostatic equilibrium.
Thermal emission and absorption obey Kirchhoff's law at the electron temperature.
Compton scattering is not included. The legacy synchrotron fit is cut off below
5e8 K; this approximation must be assessed for the intended frequency range.
"""
@with_kw mutable struct IonTorus{R, P} <: AbstractRadiativeModel
    λ::Float64 = 0.3 #Specific angular momentum dimensionless parameter
    ϵc::Float64 = 1e-17 #Central mass density in g cm^-3
    n::Float64 = 3 / 2 #Polytropic index
    Tec::Float64 = 2e9 #Central electorn temperature in Kelvin
    ξ::Float64 = 0.1 #Electron to proton temperature ratio at the center 
    β::Float64 = 0.1 #Equipartition factor
    H_abundance::Float64 = 0.75
    He_abundance::Float64 = 0.25
    rotation_sense::R = ProgradeRotation()
    radiative_process::P = Bremsstrahlung()
    self_absorption::Bool = true
    thermodynamics::Symbol = :prescribed #No geometry is available in the raw constructor.
    μᵢ::Float64 = 4 / (4H_abundance + He_abundance)
    μₑ::Float64 = 1 / (H_abundance + He_abundance / 2)
    𝓜₀::Float64 = μᵢ / (μᵢ + μₑ)
    𝓜₁::Float64 = μᵢ * ξ / (μᵢ * ξ + μₑ)
    K::Float64 = PhysicalConstants.k_B * Tec /
                 ((1 - β) * PhysicalConstants.mu * ϵc^(1 / n) * μₑ * 𝓜₁)
    Hc::Float64 = (n + 1) * log1p(K * ϵc^(1 / n) / PhysicalConstants.c^2)
    l0::Float64 = 0.0
    rcusp::Float64 = 0.0
    rcenter::Float64 = 0.0
    router::Float64 = Inf
    potential_at_surface::Float64 = 0.0
    potential_at_center::Float64 = 0.0
    _l0_is_set::Bool = false
    _radii_are_set::Bool = false
    @assert 0<λ<1 "A closed nondegenerate cusp-filled torus requires 0<lambda<1"
    @assert isfinite(ϵc) && ϵc>0 && isfinite(Tec) && Tec>0 "Density and temperature must be finite and positive"
    @assert isfinite(n) && n>0 && 0<ξ≤1 && 0≤β<1 "Require n>0, 0<xi<=1 and 0<=beta<1"
    @assert H_abundance>=0 && He_abundance>=0 && isapprox(H_abundance+He_abundance,1; atol=1e-12) "H and He mass fractions must be nonnegative and sum to one"
    @assert thermodynamics in (:hydrostatic,:prescribed) "Unknown thermodynamics mode"
    @assert isa(rotation_sense, AbstractRotationSense) "Rotation sense must be either ProgradeRotation() or RetrogradeRotation()"
end

function IonTorus(spacetime::AbstractSpacetime; thermodynamics=:hydrostatic,
    Tec=nothing, l0=nothing, kwargs...)
    spacetime isa Union{KerrSpacetimeBoyerLindquistCoordinates,FRKerrSpacetime} ||
        throw(ArgumentError("IonTorus requires Kerr Boyer-Lindquist or FRKerr coordinates; use a=0 for Schwarzschild."))
    spacetime isa FRKerrSpacetime && spacetime.R0 > 0 &&
        throw(ArgumentError("Cusp-filled de Sitter tori require an outer-barrier construction, which is not implemented."))
    model = IonTorus(; Tec=isnothing(Tec) ? 2e9 : Tec, thermodynamics, kwargs...)
    if isnothing(l0)
        torus_specific_angular_momentum!(model, spacetime)
    else
        model.l0 = l0
        model._l0_is_set = true
    end
    cusp_and_center_radius!(model, spacetime)
    torus_potentials_at_center_and_surface!(model, spacetime)
    outer_torus_radius!(model, spacetime)
    if thermodynamics == :hydrostatic
        model.Hc = model.potential_at_center-model.potential_at_surface
        model.K = PhysicalConstants.c^2 * expm1(model.Hc/(model.n+1)) / model.ϵc^(1/model.n)
        equilibrium_Tec = (1-model.β)*PhysicalConstants.mu*model.μₑ*model.𝓜₁*
            model.K*model.ϵc^(1/model.n)/PhysicalConstants.k_B
        isnothing(Tec) || isapprox(Tec,equilibrium_Tec; rtol=1e-6) ||
            throw(ArgumentError("Tec does not satisfy hydrostatic balance (expected $equilibrium_Tec K). Omit Tec or explicitly use thermodynamics=:prescribed."))
        model.Tec = equilibrium_Tec
    end
    return model
end

stationarity(::IonTorus) = IsStationary()
axisymmetry(::IonTorus) = IsAxisymmetric()

isvacuum(::IonTorus) = NonVacuum()

function torus_specific_angular_momentum!(model::IonTorus, spacetime)
    λ = model.λ
    lms = innermost_stable_specific_angular_momentum(spacetime, model.rotation_sense)
    lmb = marginally_bound_specific_angular_momentum(spacetime, model.rotation_sense)
    model.l0 = lms + λ * (lmb - lms)
    model._l0_is_set = true
    return nothing
end

function cusp_and_center_radius!(model::IonTorus, spacetime)
    if !model._l0_is_set
        throw(ArgumentError("Specific angular momentum has not been set."))
    end
    M = spacetime.M
    l0 = model.l0
    rmb = mbco_radius(spacetime, model.rotation_sense)
    rms = isco_radius(spacetime, model.rotation_sense)
    sense = sign(model.rotation_sense)
    lms = innermost_stable_specific_angular_momentum(spacetime,model.rotation_sense)
    lmb = marginally_bound_specific_angular_momentum(spacetime,model.rotation_sense)
    sense*lms < sense*l0 < sense*lmb ||
        throw(ArgumentError("l0 must lie strictly between the stable and inner E=1 reference angular momenta."))
    function aux(r)
        sense*(circular_geodesic_specific_angular_momentum([0.0, r, π / 2, 0.0],
            spacetime,
            model.rotation_sense) - l0)
    end
    model.rcusp = find_zero(aux, (rmb, rms), Bisection())
    rmax = 2rms
    while aux(rmax) <= 0 && rmax < 1e8M
        rmax *= 2
    end
    aux(rmax)>0 || throw(ArgumentError("No torus center on the outer stable branch."))
    model.rcenter = find_zero(aux, (rms, rmax), Bisection())
    model._radii_are_set = true
    return nothing
end

function constant_angular_momentum_angular_speed(g::AbstractMatrix, model::IonTorus)
    l0 = model.l0
    gtt = g[1, 1]
    gtφ = g[1, 4]
    gφφ = g[4, 4]
    return -(gtφ + l0 * gtt) / (gφφ + l0 * gtφ)
end

function torus_potential(r::Real, g::AbstractMatrix, model::IonTorus)
    r >= model.rcusp && g[4,4] > 0 || return -Inf
    gtt = g[1, 1]
    gtφ = g[1, 4]
    gφφ = g[4, 4]
    Ω = constant_angular_momentum_angular_speed(g, model)
    p2 = gtt + 2Ω * gtφ + Ω^2 * gφφ
    p2 < 0 && isfinite(p2) || return -Inf
    return 0.5 * log(-p2 / (gtt + Ω * gtφ)^2)
end

function torus_potential_at_surface(spacetime, model::IonTorus)
    r = model.rcusp
    position = equatorial_position(r, coordinates_topology(spacetime))
    g = metric(position, spacetime)
    return torus_potential(r, g, model)
end

function torus_potential_at_center(spacetime, model::IonTorus)
    r = model.rcenter
    position = equatorial_position(r, coordinates_topology(spacetime))
    g = metric(position, spacetime)
    return torus_potential(r, g, model)
end

function torus_potentials_at_center_and_surface!(model::IonTorus, spacetime)
    model.potential_at_surface = torus_potential_at_surface(spacetime, model)
    model.potential_at_center = torus_potential_at_center(spacetime, model)
    isfinite(model.potential_at_center) && isfinite(model.potential_at_surface) &&
        model.potential_at_center > model.potential_at_surface ||
        throw(ArgumentError("The torus must have a finite positive potential depth."))
    return nothing
end

function outer_torus_radius!(model::IonTorus, spacetime)
    f(r) = torus_potential([0.,r,π/2,0.],spacetime,model)-model.potential_at_surface
    rmax = 2model.rcenter
    while f(rmax) > 0 && rmax < 1e9*spacetime.M
        rmax *= 2
    end
    f(rmax) < 0 || throw(ArgumentError("The torus surface does not close."))
    model.router = find_zero(f,(model.rcenter,rmax),Bisection())
    return nothing
end

function torus_normalized_potential(r::Real, g::AbstractMatrix, model::IonTorus)
    model.rcusp <= r <= model.router || return 0.0
    W = torus_potential(r, g, model)
    return torus_normalized_potential(W, model)
end

function torus_normalized_potential(W::Real, model::IonTorus)
    Ws = model.potential_at_surface
    Wc = model.potential_at_center
    ω = clamp((W - Ws) / (Wc - Ws), 0.0, 1.0)
    return ω
end

function energy_density(ω::Real, model::IonTorus)
    ω <= 0 && return 0.0
    0 < ω <= 1 || throw(DomainError(ω,"The normalized torus potential must be in [0,1]."))
    n = model.n
    K = model.K
    ϵc = model.ϵc
    qc = K * ϵc^(1/n) / PhysicalConstants.c^2
    ϵ = ϵc * (expm1(ω*log1p(qc))/qc)^n
    return ϵ
end

function pressure(ϵ::Real, model::IonTorus)
    n = model.n
    K = model.K
    return K * ϵ^(1 + 1 / n)
end

function electron_temperature(ω::Real, ϵ::Real, model::IonTorus)
    ϵ <= 0 && return 0.0
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    K = model.K
    n = model.n
    factor = (1 - β) * mu * K * ϵ^(1 / n) / k_B
    return ((1 - ω) * 𝓜₀ + ω * 𝓜₁) * μₑ * factor
end

function ion_temperature(ω::Real, ϵ::Real, model::IonTorus)
    ϵ <= 0 && return 0.0
    mu = PhysicalConstants.mu
    k_B = PhysicalConstants.k_B
    K = model.K
    n = model.n
    β = model.β
    𝓜₀ = model.𝓜₀
    𝓜₁ = model.𝓜₁
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    factor = (1 - β) * mu * K * ϵ^(1 / n) / k_B
    return ((μₑ / μᵢ)𝓜₀ + ω * (𝓜₀ - 𝓜₁)) * μᵢ * factor
end

function number_densities(ϵ::Real, model::IonTorus)
    mu = PhysicalConstants.mu
    μₑ = model.μₑ
    μᵢ = model.μᵢ
    ne = ϵ / (mu * μₑ)
    ni = ϵ / (mu * μᵢ)
    return ne, ni
end

function electron_number_density(ϵ::Real, model::IonTorus)
    mu = PhysicalConstants.mu
    μₑ = model.μₑ
    ne = ϵ / (mu * μₑ)
    return ne
end

function magnetic_field(ϵ::Real, model::IonTorus)
    β = model.β
    P = pressure(ϵ, model)
    return sqrt(24π * β * P)
end

function number_densities_and_electron_temperature(ω::Real, model::IonTorus)
    ϵ = energy_density(ω, model)
    ne, ni = number_densities(ϵ, model)
    Te = electron_temperature(ω, ϵ, model)
    return ne, ni, Te
end

function electron_number_density_temperature_and_magnetic_field(ω::Real, model::IonTorus)
    ϵ = energy_density(ω, model)
    ne = electron_number_density(ϵ, model)
    Te = electron_temperature(ω, ϵ, model)
    B = magnetic_field(ϵ, model)
    return ne, Te, B
end

function number_densities_electron_temperature_and_magnetic_field(ω::Real, model::IonTorus)
    ϵ = energy_density(ω, model)
    ne, ni = number_densities(ϵ, model)
    Te = electron_temperature(ω, ϵ, model)
    B = magnetic_field(ϵ, model)
    return ne, ni, Te, B
end

function rest_frame_four_velocity!(vector,
    position,
    metric,
    spacetime,
    model::IonTorus,
    coords_top)
    inside = torus_normalized_potential(position[2],metric,model) > 0
    # In vacuum use a locally nonrotating observer, which stays timelike in
    # the ergosphere. Its velocity does not contribute any radiation.
    Ω = inside ? constant_angular_momentum_angular_speed(metric, model) :
        (metric[4,4] > 0 ? -metric[1,4]/metric[4,4] : 0.0)
    p2 = metric[1, 1] + 2Ω * metric[1, 4] + Ω^2 * metric[4, 4]
    if p2 < 0 && isfinite(p2)
        circular_motion_four_velocity!(vector,position,Ω,metric,coords_top)
    else
        inside && throw(DomainError(p2,"The emitting fluid must be timelike."))
        fill!(vector,0.0) #Vacuum marker for an intermediate step outside this chart.
    end
    return nothing
end

"""Charge-squared weighted ion density for fully ionized H/He, in cm^-3."""
charge_squared_ion_density(ρ, model::IonTorus) =
    (model.H_abundance+model.He_abundance)*ρ/PhysicalConstants.mu

function _torus_coefficient(::Bremsstrahlung,E,ne,nZ2,Te,B,absorption)
    return absorption ? bremsstrahlung_absorptivity(E,ne,nZ2,Te) :
        bremsstrahlung_emissivity(E,ne,nZ2,Te)
end

function _torus_coefficient(sy::Synchrotron,E,ne,nZ2,Te,B,absorption)
    jE = synchrotron_emissivity(E,ne,Te,B,sy)
    return absorption ? thermal_absorptivity(jE,E,Te) : jE
end

function _torus_coefficient(syb::SynchrotronAndBremsstrahlung,E,ne,nZ2,Te,B,absorption)
    jsy = synchrotron_emissivity(E,ne,Te,B,syb)
    sy = absorption ? thermal_absorptivity(jsy,E,Te) : jsy
    return sy + _torus_coefficient(Bremsstrahlung(),E,ne,nZ2,Te,B,absorption)
end

function _torus_coefficients!(coefficients,position,energies,g,model,absorption)
    fill!(coefficients,0.0)
    ω = torus_normalized_potential(position[2],g,model)
    ω > 0 || return nothing
    ρ = energy_density(ω,model)
    Te = electron_temperature(ω,ρ,model)
    (ρ > 0 && Te > 0) || return nothing
    ne = electron_number_density(ρ,model)
    nZ2 = charge_squared_ion_density(ρ,model)
    B = magnetic_field(ρ,model)
    for i in eachindex(energies,coefficients)
        coefficients[i] = _torus_coefficient(model.radiative_process,energies[i],ne,nZ2,Te,B,absorption)
    end
    return nothing
end

function rest_frame_absorptivity!(αε,position,ε,g,spacetime,model::IonTorus,
    coords_top,spacetime_cache,model_cache)
    if model.self_absorption
        _torus_coefficients!(αε,position,ε,g,model,true)
    else
        fill!(αε,0.0)
    end
    return nothing
end

function rest_frame_emissivity!(jε,position,ε,g,spacetime,model::IonTorus,
    coords_top,spacetime_cache,model_cache)
    _torus_coefficients!(jε,position,ε,g,model,false)
end

function energy_density(r::Real, g::AbstractMatrix, model::IonTorus)
    ω = torus_normalized_potential(r, g, model)
    return energy_density(ω, model)
end

function pressure(r::Real, g::AbstractMatrix, model::IonTorus)
    ϵ = energy_density(r, g, model)
    return pressure(ϵ, model)
end

function electron_temperature(r::Real, g::AbstractMatrix, model::IonTorus)
    ω = torus_normalized_potential(r, g, model)
    ϵ = energy_density(ω, model)
    return electron_temperature(ω, ϵ, model)
end

function ion_temperature(r::Real, g::AbstractMatrix, model::IonTorus)
    ω = torus_normalized_potential(r, g, model)
    ϵ = energy_density(ω, model)
    return ion_temperature(ω, ϵ, model)
end

function number_densities(r::Real, g::AbstractMatrix, model::IonTorus)
    ϵ = energy_density(r, g, model)
    return number_densities(ϵ, model)
end

function torus_potential(position, spacetime, model::IonTorus)
    torus_potential(position[2], metric(position, spacetime), model)
end
function torus_normalized_potential(position, spacetime, model::IonTorus)
    torus_normalized_potential(position[2], metric(position, spacetime), model)
end
function energy_density(position, spacetime, model::IonTorus)
    energy_density(position[2], metric(position, spacetime), model)
end
function pressure(position, spacetime, model::IonTorus)
    pressure(position[2], metric(position, spacetime), model)
end
function electron_temperature(position, spacetime, model::IonTorus)
    electron_temperature(position[2], metric(position, spacetime), model)
end
function ion_temperature(position, spacetime, model::IonTorus)
    ion_temperature(position[2], metric(position, spacetime), model)
end
function number_densities(position, spacetime, model::IonTorus)
    number_densities(position[2], metric(position, spacetime), model)
end
function magnetic_field(position, spacetime, model::IonTorus)
    magnetic_field(energy_density(position, spacetime, model), model)
end
function number_densities_and_electron_temperature(position, spacetime, model::IonTorus)
    number_densities_and_electron_temperature(torus_normalized_potential(position,
            spacetime,
            model),
        model)
end
function electron_number_density_temperature_and_magnetic_field(position,
    spacetime,
    model::IonTorus)
    electron_number_density_temperature_and_magnetic_field(torus_normalized_potential(position,
            spacetime,
            model),
        model)
end
