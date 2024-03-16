@doc raw"""
    FluidStarSpacetime <: AbstractSpacetime

Represents the spacetime of a fluid star, using spherical coordinates. This model includes both the internal and external spacetime metrics of a spherically symmetric, non-rotating fluid star of mass `M` and radius `R`. The internal metric (applicable for `r <= R`) accounts for the star's material properties, while the external metric (for `r > R`) follows the Schwarzschild solution, indicative of the spacetime in a vacuum outside the star.

The internal metric is described by:

\[
ds^2_- = -\frac{1}{4}\left(3\sqrt{1-\frac{2M}{R}}-\sqrt{1-\frac{2r^2M}{R^3}}\right)^2dt^2 + \left(1-\frac{2r^2M}{R^3}\right)^{-1}dr^2 + r^2\left(d\theta^2+\sin^2\theta d\phi^2\right),
\]

and the external metric by:

\[
ds^2_+ = -\left(1-\frac{2M}{r}\right)dt^2 + \left(1-\frac{2M}{r}\right)^{-1}dr^2 + r^2\left(d\theta^2+\sin^2\theta d\phi^2\right).
\]

# Fields
- `M::Float64`: The mass of the fluid star.
- `R::Float64`: The radius of the fluid star.tors
"""
@with_kw struct FluidStarSpacetime <: AbstractRegularCompactObjectSpacetime
    M::Float64  # Mass of the fluid star
    R::Float64  # Radius of the fluid star
    @assert M > 0 "The mass of the fluid star must be positive"
    @assert R > 2*M "The radius of the fluid star must be greater than 2 times its mass"
end

stationarity(::FluidStarSpacetime) = IsStationary()
spherical_symmetry(::FluidStarSpacetime) = IsSphericallySymmetric()

coordinates_topology(::FluidStarSpacetime) = SphericalTopology()
radius(position, ::FluidStarSpacetime) = position[2]

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::FluidStarSpacetime)
    r = position[2]
    θ = position[3]
    M = spacetime.M
    R = spacetime.R
    in = r <= R

    gtt_internal = -0.25 * (3 * sqrt(1 - 2 * M / R) - (sqrt∘abs)(1 - 2 * r^2 * M / R^3))^2 #abs for this to be calculated also outside, for inline if
    grr_internal = 1 / (1 - 2 * r^2 * M / R^3)
    gtt_external = -(1 - 2 * M / r)
    grr_external = 1 / (1 - 2 * M / r)
    gtt = in*gtt_internal + (1 - in)*gtt_external
    grr = in*grr_internal + (1 - in)*grr_external

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

function metric_inverse!(g::AbstractMatrix, position::AbstractVector, spacetime::FluidStarSpacetime, ::AbstractMatrix, ::Nothing)
    r = position[2]
    θ = position[3]
    M = spacetime.M
    R = spacetime.R
    in = r <= R

    gtt_internal = -0.25 * (3 * sqrt(1 - 2 * M / R) - (sqrt∘abs)(1 - 2 * r^2 * M / R^3))^2 #abs for this to be calculated also outside, for inline if
    grr_internal = 1 / (1 - 2 * r^2 * M / R^3)
    gtt_external = -(1 - 2 * M / r)
    grr_external = 1 / (1 - 2 * M / r)
    gtt = in*gtt_internal + (1 - in)*gtt_external
    grr = in*grr_internal + (1 - in)*grr_external

    g[1, 1] = 1 / gtt
    g[1, 2] = 0.0
    g[1, 3] = 0.0
    g[1, 4] = 0.0
    g[2, 1] = 0.0
    g[2, 2] = 1 / grr
    g[2, 3] = 0.0
    g[2, 4] = 0.0
    g[3, 1] = 0.0
    g[3, 2] = 0.0
    g[3, 3] = 1.0 / r^2
    g[3, 4] = 0.0
    g[4, 1] = 0.0
    g[4, 2] = 0.0
    g[4, 3] = 0.0
    g[4, 4] = 1.0 / (r^2 * sin(θ)^2)
    return nothing
end

# allocate_christoffel_cache(::FluidStarSpacetime) = nothing

# function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::FluidStarSpacetime)
#     #Spacetime coordinates
#     error("Not implemented")
#     r = position[2]
#     θ = position[3]

#     gtt = -1 + numa / dena
#     grr = 1 - numb / denb


#     ∂r_gtt = ∂r_numa / dena - numa * ∂r_dena / dena^2
#     ∂r_grr = -∂r_numb / denb + numb * ∂r_denb / denb^2

#     dα = ∂r_gtt / (2 * gtt)
#     dβ = ∂r_grr / (2 * grr)

#     # Γ[1, 1, 2] = dα
#     # Γ[1, 2, 1] = Γ[1, 1, 2]

#     # Γ[2, 1, 1] = -gtt / grr * dα
#     # Γ[2, 2, 2] = dβ
#     # Γ[2, 3, 3] = -r / grr
#     # Γ[2, 4, 4] = -r / grr * sin(θ)^2

#     # Γ[3, 2, 3] = 1.0 / r
#     # Γ[3, 4, 4] = -sin(θ) * cos(θ)
#     # Γ[3, 3, 2] = Γ[3, 2, 3]

#     # Γ[4, 2, 4] = 1.0 / r
#     # Γ[4, 3, 4] = cos(θ) / sin(θ)
#     # Γ[4, 4, 2] = Γ[4, 2, 4]
#     # Γ[4, 4, 3] = Γ[4, 3, 4]
#     return nothing
# end

function circular_geodesic_angular_speed(position,
    spacetime::FluidStarSpacetime,
    rotation_sense)
    #Spacetime coordinates
    r = position[2]
    M = spacetime.M
    Ω = sqrt(M) / r^1.5
    s = sign(rotation_sense)
    return s * Ω
end