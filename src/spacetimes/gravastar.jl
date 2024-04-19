@doc raw"""
    GravastarSpacetime <: AbstractRegularCompactObjectSpacetime

Represents the spacetime of a gravastar, using spherical coordinates. This model includes both the internal and external spacetime metrics of a spherically symmetric gravastar. The internal metric is characterized by a constant \(\alpha\) that influences the gravastar's effective mass distribution, and \(M_\rho\), the volumetric mass. The external metric follows the Schwarzschild solution, indicative of the spacetime in a vacuum outside the gravastar.

The internal metric is described by:

\[
ds^2_- = -\alpha\left(1-\frac{2r^2M_\rho}{R^3}\right)dt^2 + \left(1-\frac{2r^2M_\rho}{R^3}\right)^{-1}dr^2 + r^2\left(d\theta^2+\sin^2\theta d\phi^2\right),
\]

and the external metric by:

\[
ds^2_+ = -\left(1-\frac{2M}{r}\right)dt^2 + \left(1-\frac{2M}{r}\right)^{-1}dr^2 + r^2\left(d\theta^2+\sin^2\theta d\phi^2\right).
\]

The parameter \(\alpha\) is defined as:

\[
\alpha = \frac{1-\frac{2M}{R}}{1-\frac{2M_\rho}{R}},
\]

indicating the distribution of mass within the gravastar. This model is described by the two free parameters \(\alpha\) and \(R\), with constraints to ensure physical relevance.

# Fields
- `M::Float64`: The mass of the gravastar.
- `M_rho::Float64`: The volumetric mass of the gravastar.
- `R::Float64`: The radius of the gravastar.
- `alpha::Float64`: The parameter controlling the mass distribution.
"""
@with_kw struct GravastarSpacetime <: AbstractRegularCompactObjectSpacetime
    M::Float64
    M_rho::Float64
    R::Float64
    @assert M > 0 "The mass of the gravastar must be positive"
    @assert M_rho > 0 "The volumetric mass of the gravastar must be positive"
    @assert R > 2*M "The radius of the gravastar must be greater than 2 times its mass"
    alpha::Float64 = (1 - 2*M / R) / (1 - 2*M_rho / R)
    @assert alpha >= (1 - 2*M/R) "Alpha must not lead to negative energy densities, ensuring physical relevance"
end
    
stationarity(::GravastarSpacetime) = IsStationary()
spherical_symmetry(::GravastarSpacetime) = IsSphericallySymmetric()

coordinates_topology(::GravastarSpacetime) = SphericalTopology()
radius(position, ::GravastarSpacetime) = position[2]

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::GravastarSpacetime)
    r = position[2]
    θ = position[3]
    M = spacetime.M
    M_rho = spacetime.M_rho
    R = spacetime.R
    alpha = spacetime.alpha
    in = r <= R

    gtt_internal = -alpha * (1 - 2 * r^2 * M_rho / R^3)
    gtt_external = -(1 - 2 * M / r)

    gtt = in * gtt_internal + (1 - in) * gtt_external
    grr = -1/(gtt*(in/alpha + (1 - in)))

    # The spherical symmetry parts remain unchanged
    g[1, 1] = gtt
    g[1, 2] = 0.0
    g[1, 3] = 0.0
    g[1, 4] = 0.0
    g[2, 1] = 0.0
    g[2, 2] = grr
    g[2, 3] = 0.0
    g[2, 4] = 0.0
    g[3, 1] = 0.0
    g[3, 2] = 0.0
    g[3, 3] = r^2
    g[3, 4] = 0.0
    g[4, 1] = 0.0
    g[4, 2] = 0.0
    g[4, 3] = 0.0
    g[4, 4] = r^2 * sin(θ)^2
    return nothing
end
    
allocate_christoffel_cache(::GravastarSpacetime) = nothing

function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::GravastarSpacetime)
    #Spacetime coordinates
    r = position[2]
    θ = position[3]
    M = spacetime.M
    M_rho = spacetime.M_rho
    R = spacetime.R
    alpha = spacetime.alpha
    in = r <= R

    gtt_internal = -alpha * (1 - 2 * r^2 * M_rho / R^3)
    gtt_external = -(1 - 2 * M / r)

    gtt = in * gtt_internal + (1 - in) * gtt_external
    grr = -1/(gtt*(in/alpha + (1 - in)))

    ∂r_gtt_internal = alpha*4*r*M_rho/R^3
    ∂r_gtt_external = -2 * M / r^2

    ∂r_gtt = in*∂r_gtt_internal + (1 - in)*∂r_gtt_external
    ∂r_grr = (-∂r_gtt / gtt^2)*(in*alpha + (1 - in))
    
    dα = ∂r_gtt / (2 * gtt)
    dβ = ∂r_grr / (2 * grr)

    Γ[1, 1, 2] = dα
    Γ[1, 2, 1] = Γ[1, 1, 2]

    Γ[2, 1, 1] = -gtt / grr * dα
    Γ[2, 2, 2] = dβ
    Γ[2, 3, 3] = -r / grr
    Γ[2, 4, 4] = -r / grr * sin(θ)^2

    Γ[3, 2, 3] = 1.0 / r
    Γ[3, 4, 4] = -sin(θ) * cos(θ)
    Γ[3, 3, 2] = Γ[3, 2, 3]

    Γ[4, 2, 4] = 1.0 / r
    Γ[4, 3, 4] = cos(θ) / sin(θ)
    Γ[4, 4, 2] = Γ[4, 2, 4]
    Γ[4, 4, 3] = Γ[4, 3, 4]
    return nothing
end

function circular_geodesic_angular_speed(position,
    spacetime::GravastarSpacetime,
    rotation_sense)
    M = spacetime.M
    r = radius(position, spacetime)
    Ω = sqrt(M) / r^1.5
    s = sign(rotation_sense)
    return s * Ω
end

mass(spacetime::GravastarSpacetime) = spacetime.M

function isco_radius(spacetime::GravastarSpacetime, ::AbstractRotationSense) 
    spacetime.R <= 6*spacetime.M || error(ArgumentError("The radius of the star is larger than 6M. The ISCO is not implemented."))
    return 6 * spacetime.M
end