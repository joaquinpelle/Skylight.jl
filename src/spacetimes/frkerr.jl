@doc raw"""
    FRKerrSpacetime <: AbstractBlackHoleSpacetime

[$f(R)$-Kerr spacetime](https://www.worldscientific.com/doi/abs/10.1142/S0219887814500017) 
(uncharged) with constant curvature $R_0$

# Constructor
```
FRKerrSpacetime(M=1.0, a=0.99, R0=-0.0012) 
```
"""
@with_kw struct FRKerrSpacetime <: AbstractBlackHoleSpacetime
    M::Float64
    a::Float64
    R0::Float64
    @assert isfinite(M) && M>0.0 "M must be finite and positive"
    @assert isfinite(a) && isfinite(R0) "a and R0 must be finite"
    @assert 1 + a^2*R0/12 > 0 "Xi must be positive"
end

stationarity(::FRKerrSpacetime) = IsStationary()
axisymmetry(::FRKerrSpacetime) = IsAxisymmetric()

coordinates_topology(::FRKerrSpacetime) = SphericalTopology()
radius(position, ::FRKerrSpacetime) = position[2]

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::FRKerrSpacetime)
    r = position[2]
    θ = position[3]
    M = spacetime.M
    a = spacetime.a
    R0 = spacetime.R0

    ρ² = r^2 + a^2 * cos(θ)^2
    Δᵣ = (r^2 + a^2) * (1 - R0 * r^2 / 12) - 2 * M * r
    Δθ = 1 + R0 * a^2 * cos(θ)^2 / 12
    Ξ = 1 + a^2 * R0 / 12

    fill!(g, 0.0)
    g[1, 1] = (a^2 * Δθ * sin(θ)^2 - Δᵣ) / (ρ² * Ξ^2)
    g[1, 4] = (-a * (r^2 + a^2) * Δθ + a * Δᵣ) * sin(θ)^2 / (ρ² * Ξ^2)
    g[2, 2] = ρ² / Δᵣ
    g[3, 3] = ρ² / Δθ
    g[4, 1] = g[1, 4]
    g[4, 4] = ((r^2 + a^2)^2 * Δθ - a^2 * sin(θ)^2 * Δᵣ) * sin(θ)^2 / (ρ² * Ξ^2)
    return nothing
end

# Work with the original polynomial, which also has a regular Kerr limit.
_fr_delta(r, s::FRKerrSpacetime) = (r^2 + s.a^2) * (1 - s.R0*r^2/12) - 2*s.M*r
_fr_delta_prime(r, s::FRKerrSpacetime) =
    2r * (1 - s.R0*r^2/12) - s.R0*r*(r^2+s.a^2)/6 - 2*s.M

function horizon_parameter(s::FRKerrSpacetime)
    iszero(s.R0) && throw(ArgumentError("The horizon discriminant parameter is undefined at R0=0; use horizons()."))
    a, R0, M = s.a, s.R0, s.M
    return (4/R0*(1-R0*a^2/12)^2-4a^2)^3 +
        4/R0*((1-R0*a^2/12)*(4/R0*(1-R0*a^2/12)^2+12a^2)-18M^2)^2
end

"""Sorted real roots of Δr, including negative roots and repeated Kerr roots."""
function horizons(s::FRKerrSpacetime)
    if iszero(s.R0)
        d = s.M^2 - s.a^2
        return d < 0 ? Float64[] : [s.M-sqrt(d), s.M+sqrt(d)]
    end
    # Dimensionless radii improve conditioning when M is not one.
    a, y = s.a/s.M, s.R0*s.M^2/12
    xi = roots(Polynomial([a^2, -2.0, 1-y*a^2, 0.0, -y]))
    tol = 32sqrt(eps(Float64))
    return sort!([s.M*real(x) for x in xi if abs(imag(x)) <= tol*max(1,abs(real(x)))])
end

mass(s::FRKerrSpacetime) = s.M
spin(s::FRKerrSpacetime) = s.a

function event_horizon_radius(s::FRKerrSpacetime)
    if iszero(s.R0)
        abs(s.a) <= s.M || throw(ArgumentError("No black-hole horizon."))
        return s.M + sqrt(s.M^2 - s.a^2)
    end
    # The event horizon is the outward crossing into Δr>0. A lone de Sitter
    # root is cosmological and must never be used as an inner ray cutoff.
    candidates = filter(r -> r > 0 && _fr_delta_prime(r,s) > 1e-7*s.M, horizons(s))
    length(candidates) == 1 || throw(ArgumentError("No unique nondegenerate black-hole event horizon; inspect horizons()."))
    return only(candidates)
end

function cosmologic_horizon_radius(s::FRKerrSpacetime)
    s.R0 > 0 || throw(ArgumentError("A cosmological horizon requires R0>0."))
    rh = event_horizon_radius(s)
    candidates = filter(r -> r > rh && _fr_delta_prime(r,s) < -1e-7*s.M, horizons(s))
    length(candidates) == 1 || throw(ArgumentError("No unique nondegenerate cosmological horizon."))
    return only(candidates)
end

function _fr_orbit_terms(r, spacetime, rotation_sense)
    x = r / spacetime.M
    a = spacetime.a / spacetime.M
    y = spacetime.R0 * spacetime.M^2 / 12
    s = sign(rotation_sense)
    z = sqrt(x^(-3) - y)
    return x, a, y, s, z
end

"""Rescaled circular-orbit energy E=-Ξu_t (not -u_t when R0≠0)."""
function circular_geodesic_energy(position, spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    x,a,y,s,z = _fr_orbit_terms(radius(position,spacetime),spacetime,rotation_sense)
    return (1-2/x-(x^2+a^2)*y+s*a*z)/sqrt(1-3/x-a^2*y+2s*a*z)
end

function circular_geodesic_specific_angular_momentum(position, spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    x,a,y,s,z = _fr_orbit_terms(radius(position,spacetime),spacetime,rotation_sense)
    return spacetime.M * (s*x*(x^2+a^2)*z-2a-a*x*(x^2+a^2)*y) /
        (x*(1-2/x-(x^2+a^2)*y+s*a*z))
end

function circular_geodesic_angular_speed(position, spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    x,a,y,s,z = _fr_orbit_terms(radius(position,spacetime),spacetime,rotation_sense)
    return s*z / (spacetime.M*(1+s*a*z))
end

function _fr_bracket(f, lo, hi; rising=true)
    hi > lo || throw(ArgumentError("Empty physical circular-orbit domain."))
    previous_x, previous_y = lo, f(lo)
    for x in exp.(range(log(lo), log(hi); length=2049))[2:end]
        y = f(x)
        crossing = rising ? previous_y <= 0 <= y : previous_y >= 0 >= y
        if isfinite(previous_y) && isfinite(y) && crossing
            return previous_x, x
        end
        previous_x, previous_y = x, y
    end
    throw(ArgumentError("The requested orbit does not exist on the physical exterior branch."))
end

"""Outer reality limit of circular geodesics; this is not the OSCO."""
function outermost_circular_orbit_radius(s::FRKerrSpacetime, ::AbstractRotationSense)
    return s.R0 > 0 ? cbrt(12*s.M/s.R0) : Inf
end

function _fr_search_limit(s, sense)
    return s.R0 > 0 ? min(outermost_circular_orbit_radius(s,sense),
        cosmologic_horizon_radius(s))*(1-1e-8) : 1e3*s.M
end

function innermost_circular_orbit_radius(spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    rh = event_horizon_radius(spacetime)
    f(r) = begin
        x,a,y,s,z = _fr_orbit_terms(r,spacetime,rotation_sense)
        1-3/x-a^2*y+2s*a*z
    end
    bracket = _fr_bracket(f,rh*(1+1e-8),_fr_search_limit(spacetime,rotation_sense))
    return find_zero(f,bracket,Bisection())
end

function isco_radius(spacetime::FRKerrSpacetime, rotation_sense::AbstractRotationSense)
    rph = innermost_circular_orbit_radius(spacetime,rotation_sense)
    # For the negative-angular-momentum branch, minimize its magnitude, not l.
    f(r) = ForwardDiff.derivative(r -> sign(rotation_sense) *
        circular_geodesic_specific_angular_momentum([0.,r,π/2,0.],spacetime,rotation_sense),r)
    bracket = _fr_bracket(f,rph*(1+1e-7),_fr_search_limit(spacetime,rotation_sense))
    return find_zero(f,bracket,Bisection())
end

"""
Inner E=1 orbit, continued from Kerr. For R0<0 this is a parameterization
convention, not an escape threshold at infinity. For R0>0 the marginally bound
torus requires equal inner/outer potential barriers; that construction is not
implemented and is explicitly rejected instead of returning a spurious E=1 root.
"""
function mbco_radius(spacetime::FRKerrSpacetime, rotation_sense::AbstractRotationSense)
    spacetime.R0 <= 0 || throw(ArgumentError("For R0>0, E=1 is not the marginally bound condition; the two-barrier construction is not implemented."))
    rph = innermost_circular_orbit_radius(spacetime,rotation_sense)
    rms = isco_radius(spacetime,rotation_sense)
    f(r) = circular_geodesic_energy([0.,r,π/2,0.],spacetime,rotation_sense)-1
    lo = rph*(1+1e-7)
    f(lo) > 0 && f(rms) < 0 || throw(ArgumentError("No inner E=1 orbit outside the photon orbit and inside the ISCO."))
    return find_zero(f,(lo,rms),Bisection())
end

function innermost_stable_specific_angular_momentum(s::FRKerrSpacetime,sense)
    return circular_geodesic_specific_angular_momentum([0.,isco_radius(s,sense),π/2,0.],s,sense)
end

function marginally_bound_specific_angular_momentum(s::FRKerrSpacetime,sense)
    return circular_geodesic_specific_angular_momentum([0.,mbco_radius(s,sense),π/2,0.],s,sense)
end
