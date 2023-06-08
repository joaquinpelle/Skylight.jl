@with_kw struct FRKerrSpacetime <: AbstractBlackHoleSpacetime
    M::Float64
    a::Float64
    R0::Float64

    @assert M >= 0.0 "M must be non-negative"
end

coordinates_topology(::FRKerrSpacetime) = SphericalTopology()
radius(position, ::FRKerrSpacetime) = position[2]

function metric!(g, position, spacetime::FRKerrSpacetime)
    t, r, θ, φ = position
    M = spacetime.M
    a = spacetime.a
    R0 = spacetime.R0

    ρ² = r^2 + a^2*cos(θ)^2
    Δᵣ = (r^2 + a^2)*(1 - R0*r^2/12) - 2*M*r
    Δθ = 1 + R0*a^2*cos(θ)^2/12
    Ξ = 1 + a^2*R0/12
    
    fill!(g, 0.0)
    g[1,1] = (a^2*Δθ*sin(θ)^2-Δᵣ)/(ρ²*Ξ^2)
    g[1,4] = (-a*(r^2+a^2)*Δθ+a*Δᵣ)*sin(θ)^2/(ρ²*Ξ^2)
    g[2,2] = ρ²/Δᵣ
    g[3,3] = ρ²/Δθ
    g[4,1] = g[1,4] 
    g[4,4] = ((r^2+a^2)^2*Δθ-a*Δᵣ)*sin(θ)^2/(ρ²*Ξ^2)
    return nothing 
end

function horizon_parameter(spacetime::FRKerrSpacetime)
    M = spacetime.M
    a = spacetime.a
    R0 = spacetime.R0

    term1 = (4/R0 * (1 - R0/12 * a^2)^2 - 4 * a^2)^3
    term2 = 4/R0 * ((1 - R0/12 * a^2) * (4/R0 * (1 - R0/12 * a^2)^2 + 12 *a^2) - 18 * M^2)^2
    h = term1 + term2
    return h
end

"""May return negative roots"""
function horizons(spacetime::FRKerrSpacetime)

    R0 = spacetime.R0

    M = spacetime.M
    a = spacetime.a
    R0 = spacetime.R0

    c0 = -12/R0*a^2
    c1 = 24M/R0
    c2 = a^2-12/R0

    P = Polynomial([c0,c1,c2,0.0,1.0])

    xi = roots(P)
    rxi = xi[isreal.(xi)]
    return Real.(rxi)
end

function event_horizon_radius(spacetime::FRKerrSpacetime)

    h = horizon_parameter(spacetime)
    
    if h==0.0
        error("Horizons are degenerate. Please use horizons() to find the horizons and decide.")
    end

    ri = horizons(spacetime)
    ri = ri[ri.>0.0]

    if length(ri)==1
        return ri[1]
    elseif (length(ri)==2 || length(ri)==3)
        return ri[2]
    elseif length(ri)==0
        error("No horizons found. Please check the horizon_parameter() and horizons() functions.")
    else
        error("Something went wrong. Please check the horizon_parameter() and horizons() functions.")
    end
end

function cosmologic_horizon_radius(spacetime::FRKerrSpacetime)

    R0 = spacetime.R0
    h = horizon_parameter(spacetime)
    
    if h==0.0
        error("Horizons are degenerate. Please use horizons() to find the horizons and decide.")
    end

    ri = horizons(spacetime)
    ri = ri[ri.>0.0]

    if length(ri)==3
        return ri[3]
    else
        error("No cosmologic hoirzon. Please check the horizon_parameter() and horizons() functions.")
    end
end

"""Provisory taken equal to Kerr"""
function isco_radius(spacetime::FRKerrSpacetime, rotation_sense::AbstractRotationSense)
    x = spacetime.a/spacetime.M      #Rotation parameter
    Z1=1+cbrt(1-x^2)*(cbrt(1+x)+cbrt(1-x))
    Z2=sqrt(3*x^2+Z1^2)
    s = sign(rotation_sense)
    return spacetime.M*(3+Z2 - s*sqrt((3-Z1)*(3+Z1+2*Z2)))
end

"""Provisory taken equal to Kerr"""
function circular_geodesic_angular_speed(position, spacetime::FRKerrSpacetime, rotation_sense::AbstractRotationSense)
    M = spacetime.M
    a = spacetime.a
    r = radius(position, spacetime)
    s = sign(rotation_sense)
    Ω = s*sqrt(M)/(r^1.5 + s*a*sqrt(M))
    return Ω
end