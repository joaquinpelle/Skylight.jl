@with_kw struct FRKerrSpacetime <: AbstractBlackHoleSpacetime
    M::Float64
    a::Float64
    R0::Float64
    @assert M>=0.0 "M must be non-negative"
end

coordinates_topology(::FRKerrSpacetime) = SphericalTopology()
radius(position, ::FRKerrSpacetime) = position[2]

function metric!(g, position, spacetime::FRKerrSpacetime)
    t, r, θ, φ = position
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
    g[4, 4] = ((r^2 + a^2)^2 * Δθ - a * Δᵣ) * sin(θ)^2 / (ρ² * Ξ^2)
    return nothing
end

function horizon_parameter(spacetime::FRKerrSpacetime)
    M = spacetime.M
    a = spacetime.a
    R0 = spacetime.R0

    term1 = (4 / R0 * (1 - R0 / 12 * a^2)^2 - 4 * a^2)^3
    term2 = 4 / R0 *
            ((1 - R0 / 12 * a^2) * (4 / R0 * (1 - R0 / 12 * a^2)^2 + 12 * a^2) - 18 * M^2)^2
    h = term1 + term2
    return h
end

"""May return negative roots"""
function horizons(spacetime::FRKerrSpacetime)
    R0 = spacetime.R0

    M = spacetime.M
    a = spacetime.a
    R0 = spacetime.R0

    c0 = -12 / R0 * a^2
    c1 = 24M / R0
    c2 = a^2 - 12 / R0

    P = Polynomial([c0, c1, c2, 0.0, 1.0])

    xi = roots(P)
    rxi = xi[isreal.(xi)]
    return Real.(rxi)
end

function event_horizon_radius(spacetime::FRKerrSpacetime)
    #The horizon parameter is undefined for R0=0, so we use Kerr in that case
    if spacetime.R0 == 0.0
        return event_horizon_radius(KerrSpacetimeBoyerLindquistCoordinates(M = spacetime.M,
            a = spacetime.a))
    end

    h = horizon_parameter(spacetime)
    if h == 0.0
        error("Horizons are degenerate. Please use horizons() to find the horizons and decide.")
    end
    ri = horizons(spacetime)
    ri = ri[ri .> 0.0]
    if length(ri) == 1
        return ri[1]
    elseif (length(ri) == 2 || length(ri) == 3)
        return ri[2]
    elseif length(ri) == 0
        error("No horizons found. Please check the horizon_parameter() and horizons() functions.")
    else
        error("Something went wrong. Please check the horizon_parameter() and horizons() functions.")
    end
end

function cosmologic_horizon_radius(spacetime::FRKerrSpacetime)
    h = horizon_parameter(spacetime)
    if h == 0.0
        error("Horizons are degenerate. Please use horizons() to find the horizons and decide.")
    end
    ri = horizons(spacetime)
    ri = ri[ri .> 0.0]
    if length(ri) == 3
        return ri[3]
    else
        error("No cosmologic hoirzon. Please check the horizon_parameter() and horizons() functions.")
    end
end

"""Set for M=1"""
function circular_geodesic_energy(position,
    spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    a = spacetime.a
    R0 = spacetime.R0
    r = radius(position, spacetime)
    s = sign(rotation_sense)
    term1 = 1 - 2 / r - (r^2 + a^2) * (R0 / 12) + s * a * (1 / r^3 - (R0 / 12))^(1 / 2)
    term2 = (1 - 3 / r - a^2 * (R0 / 12) +
             s * 2 * a * (1 / r^3 - (R0 / 12))^(1 / 2))^(1 / 2)
    E = term1 / term2
    return E
end

"""Set for M=1"""
function circular_geodesic_specific_angular_momentum(position,
    spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    a = spacetime.a
    R0 = spacetime.R0
    y = R0 / 12
    r = radius(position, spacetime)
    s = sign(rotation_sense)
    A = 2 * a + a * r * (r^2 + a^2) * y - s * r * (r^2 + a^2) * (r^(-3) - y)^(1 / 2)
    B = r * (1 - 2 / r - (r^2 + a^2) * y + s * a * (r^(-3) - y)^(1 / 2))
    l_K = -A / B
    return l_K
end

"""Set for M=1"""
function isco_radius(spacetime::FRKerrSpacetime, rotation_sense::AbstractRotationSense)
    bounds = (1.0, 10.0)  # Search interval between 1 and 10 for r
    result = optimize(r -> circular_geodesic_specific_angular_momentum([0.0, r, π / 2, 0.0],
            spacetime,
            rotation_sense),
        bounds[1],
        bounds[2])
    r_isco = Optim.minimizer(result)
    return r_isco
end

"""Set for M=1"""
function mbco_radius(spacetime::FRKerrSpacetime, rotation_sense::AbstractRotationSense)
    rmin = innermost_circular_orbit_radius(spacetime, rotation_sense)
    function equation(r)
        circular_geodesic_energy([0.0, r, π / 2, 0.0], spacetime, rotation_sense) - 1
    end
    roots = find_zeros(equation, rmin + 1e-8, 10)
    r_mbco = maximum(roots)
    return r_mbco
end

function innermost_circular_orbit_radius(spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    a = spacetime.a
    R0 = spacetime.R0
    s = sign(rotation_sense)
    f(r) = (1 - 3 / r - a^2 * (R0 / 12) + s * 2 * a * (1 / r^3 - (R0 / 12))^(1 / 2))
    roots = find_zeros(f, 1, 10)
    rmin = minimum(roots)
    return rmin
end

"""Just return 100M, this if for ion torus"""
function outermost_circular_orbit_radius(spacetime::FRKerrSpacetime,
    rotation_sense::AbstractRotationSense)
    if spacetime.R0 <= 0.0
        return 10 * spacetime.M
    else
        return cbrt(12 / spacetime.R0)
    end

    return
end

function innermost_stable_specific_angular_momentum(spacetime::FRKerrSpacetime,
    rotation_sense)
    risco = isco_radius(spacetime, rotation_sense)
    position = equatorial_position(risco, coordinates_topology(spacetime))
    return circular_geodesic_specific_angular_momentum(position, spacetime, rotation_sense)
end

function marginally_bound_specific_angular_momentum(spacetime::FRKerrSpacetime,
    rotation_sense)
    rmbco = mbco_radius(spacetime, rotation_sense)
    position = equatorial_position(rmbco, coordinates_topology(spacetime))
    return circular_geodesic_specific_angular_momentum(position, spacetime, rotation_sense)
end
