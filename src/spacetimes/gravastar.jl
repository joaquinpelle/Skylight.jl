using Parameters  # For the @with_kw macro

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
        grr = 1/gtt
    
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