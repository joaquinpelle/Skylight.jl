@doc raw"""
    AxionicBosonStarSpacetime <: AbstractRegularCompactObjectSpacetime

 [Axionic boson star spacetime]() in spherical coordinates. It uses analyical fits. Either the
 fit parameters can be provided as vectors, or any of the symbols `:ABS6`, `:ABS7`, `:ABS8` as constructor arguments

# Constructors
```
AxionicBosonStarSpacetime(a=aparams, b=bparams, μ=μ) 
AxionicBosonStarSpacetime(:ABS7)
```
"""
@with_kw struct AxionicBosonStarSpacetime <: AbstractRegularCompactObjectSpacetime
    a::Vector{Float64}
    b::Vector{Float64}
    μ::Float64 

    @assert length(a)==14 "a must be a vector of length 14"
    @assert length(b)==14 "b must be a vector of length 14"
end

stationarity(::AxionicBosonStarSpacetime) = IsStationary()
spherical_symmetry(::AxionicBosonStarSpacetime) = IsSphericallySymmetric()

coordinates_topology(::AxionicBosonStarSpacetime) = SphericalTopology()
radius(position, ::AxionicBosonStarSpacetime) = position[2]

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::AxionicBosonStarSpacetime)
    r = position[2]
    θ = position[3]

    a = spacetime.a
    b = spacetime.b
    μ = spacetime.μ

    rr = r*μ

    gtt = -1 +
          (1 + rr * a[1] + rr^2 * a[2] + rr^3 * a[3] + rr^4 * a[4] + rr^5 * a[5] + rr^6 * a[6]) /
          (a[7] + rr * a[8] + rr^2 * a[9] + rr^3 * a[10] + rr^4 * a[11] + rr^5 * a[12] +
           rr^6 * a[13] + rr^7 * a[14])
    grr = 1 -
          (rr * b[1] + rr^2 * b[2] + rr^3 * b[3] + rr^4 * b[4] + rr^5 * b[5] + rr^6 * b[6]) /
          (b[7] + rr * b[8] + rr^2 * b[9] + rr^3 * b[10] + rr^4 * b[11] + rr^5 * b[12] +
           rr^6 * b[13] + rr^7 * b[14])

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

function metric_inverse!(g::AbstractMatrix, position::AbstractVector, spacetime::AxionicBosonStarSpacetime, gaux::AbstractMatrix, cache::Nothing)
    r = position[2]
    θ = position[3]

    a = spacetime.a
    b = spacetime.b
    μ = spacetime.μ
    
    rr = r*μ

    gtt = -1 +
          (1 + rr * a[1] + rr^2 * a[2] + rr^3 * a[3] + rr^4 * a[4] + rr^5 * a[5] + rr^6 * a[6]) /
          (a[7] + rr * a[8] + rr^2 * a[9] + rr^3 * a[10] + rr^4 * a[11] + rr^5 * a[12] +
           rr^6 * a[13] + rr^7 * a[14])
    grr = 1 -
          (rr * b[1] + rr^2 * b[2] + rr^3 * b[3] + rr^4 * b[4] + rr^5 * b[5] + rr^6 * b[6]) /
          (b[7] + rr * b[8] + rr^2 * b[9] + rr^3 * b[10] + rr^4 * b[11] + rr^5 * b[12] +
           rr^6 * b[13] + rr^7 * b[14])

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

allocate_christoffel_cache(::AxionicBosonStarSpacetime) = nothing

function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::AxionicBosonStarSpacetime)
    r = position[2]
    θ = position[3]

    a = spacetime.a
    b = spacetime.b
    μ = spacetime.μ

    rr = r*μ

    numa = 1 + rr * a[1] + rr^2 * a[2] + rr^3 * a[3] + rr^4 * a[4] + rr^5 * a[5] + rr^6 * a[6]
    dena = a[7] + rr * a[8] + rr^2 * a[9] + rr^3 * a[10] + rr^4 * a[11] + rr^5 * a[12] +
           rr^6 * a[13] + rr^7 * a[14]

    numb = rr * b[1] + rr^2 * b[2] + rr^3 * b[3] + rr^4 * b[4] + rr^5 * b[5] + rr^6 * b[6]
    denb = b[7] + rr * b[8] + rr^2 * b[9] + rr^3 * b[10] + rr^4 * b[11] + rr^5 * b[12] +
           rr^6 * b[13] + rr^7 * b[14]

    gtt = -1 + numa / dena
    grr = 1 - numb / denb

    ∂rr_numa = a[1] + 2rr * a[2] + 3rr^2 * a[3] + 4rr^3 * a[4] + 5rr^4 * a[5] + 6rr^5 * a[6]
    ∂rr_dena = a[8] + 2rr * a[9] + 3rr^2 * a[10] + 4rr^3 * a[11] + 5rr^4 * a[12] + 6rr^5 * a[13] +
              7rr^6 * a[14]

    ∂rr_numb = b[1] + 2rr * b[2] + 3rr^2 * b[3] + 4rr^3 * b[4] + 5rr^4 * b[5] + 6rr^5 * b[6]
    ∂rr_denb = b[8] + 2rr * b[9] + 3rr^2 * b[10] + 4rr^3 * b[11] + 5rr^4 * b[12] + 6rr^5 * b[13] +
              7rr^6 * b[14]

    ∂r_numa = μ * ∂rr_numa
    ∂r_dena = μ * ∂rr_dena
    ∂r_numb = μ * ∂rr_numb
    ∂r_denb = μ * ∂rr_denb
    
    ∂r_gtt = ∂r_numa / dena - numa * ∂r_dena / dena^2
    ∂r_grr = -∂r_numb / denb + numb * ∂r_denb / denb^2

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
    spacetime::AxionicBosonStarSpacetime,
    rotation_sense)
    r = position[2]

    a = spacetime.a
    μ = spacetime.μ
    
    rr = r*μ

    numa = 1 + rr * a[1] + rr^2 * a[2] + rr^3 * a[3] + rr^4 * a[4] + rr^5 * a[5] + rr^6 * a[6]
    dena = a[7] + rr * a[8] + rr^2 * a[9] + rr^3 * a[10] + rr^4 * a[11] + rr^5 * a[12] +
           rr^6 * a[13] + rr^7 * a[14]

    ∂rr_numa = a[1] + 2rr * a[2] + 3rr^2 * a[3] + 4rr^3 * a[4] + 5rr^4 * a[5] + 6rr^5 * a[6]
    ∂rr_dena = a[8] + 2rr * a[9] + 3rr^2 * a[10] + 4rr^3 * a[11] + 5rr^4 * a[12] + 6rr^5 * a[13] +
                7rr^6 * a[14]

    ∂r_numa = μ * ∂rr_numa
    ∂r_dena = μ * ∂rr_dena

    ∂r_gtt = ∂r_numa / dena - numa * ∂r_dena / dena^2
    ∂r_gφφ = 2r

    Ω = sqrt(-∂r_gtt / ∂r_gφφ)

    return sign(rotation_sense) * Ω
end

function AxionicBosonStarSpacetime(name::Symbol)
    
    if name==:ABS6
            a = [ -0.6494497681749755,  0.45957735938259053, 
    -0.034900123885242273,  0.10334243747951968, 
    -0.05480171014169698,  0.05442633741074079, 
    1.3619989556162608,  -0.8836727376087928, 
    0.5575020430741733,  0.04067746568545405, 
    0.028958171701746225,  0.08153232458605991, 
    -0.04584970391123078,  0.04578614418713543]

        b = [ 0.10383430635944844,  -38.56917116020277, 
        47.45976495509472,  -29.719162047716097, 
        7.13655147960732,  -1.2403548579590058, 
        73.44753758416398,  -77.61806743154906, 
        -4.554629577808198,  79.35821745158775, 
        -71.14974957902932,  32.3824545363186, 
        -7.258569429973404,  1.0438265816510337]

        μ = 0.5942940110388295

    elseif name==:ABS7

        a = [ -4.776868426909388,  10.222569885887948, 
 -15.870146868834242,  25.46143589413741, 
 -32.669611532380564,  19.4164502929398, 
 1.2169346655131477,  -5.748524999671183, 
 10.51518521245094,  -7.137076798426266, 
 -8.55167936394653,  34.47105861169818, 
 -51.614773548831415,  31.36265429822439]

b = [ -0.024144732671384544,  -2.9933921983096465, 
 20.07588687151928,  -51.85871747268101, 
 61.90268668278332,  -28.83334465925838, 
 1.1457754390531736,  -8.354912123901801, 
 22.228975979148803,  -15.726524647855072, 
 -46.9878570174763,  126.6138914354964, 
 -124.03814495147577,  46.21269047087403]

    μ = 0.3086490200937825

    elseif name==:ABS8

        a = [ -8.192553626698524,  30.601007122123235, 
 -87.04775269740725,  251.75483828432792, 
 -543.3781520837146,  526.1725491036841, 
 1.2089609070954486,  -9.731922408721891, 
 29.985241865720557,  -31.210239022301277, 
 -97.53685772912009,  589.8253269174509, 
 -1436.4312769224523,  1414.245191300103]


b = [ -0.0013750710716961322,  -0.23890122904958874, 
 2.730238660663305,  -11.856631375079349, 
 23.58360278604608,  -18.160436581818107, 
 0.03283942542232639,  -0.41224415979173484, 
 1.8992587041089353,  -2.5444475712483023, 
 -9.855263291655318,  47.21219099424397, 
 -77.93916162757478,  48.314042509721396]

    μ = 0.1946586486466792

    else
        error("Unknown name $name")
    end
    return AxionicBosonStarSpacetime(a, b, μ)
end
