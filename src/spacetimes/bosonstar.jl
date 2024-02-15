@doc raw"""
    BosonStarSpacetime <: AbstractRegularCompactObjectSpacetime

 [Boson star spacetime]() in spherical coordinates. It uses analyical fits. Either the
 fit parameters as described in... can be provided as vectors, or any of the symbols `:LBS1`,
 `:LBS2`, `:LBS3`, `:SBS1`, `:SBS2` or `:SBS3` as constructor arguments

# Constructors
```
BosonStarSpacetime(a=aparams, b=bparams) 
BosonStarSpacetime(:LBS1)
```
"""
@with_kw struct BosonStarSpacetime <: AbstractRegularCompactObjectSpacetime
    a::Vector{Float64}
    b::Vector{Float64}

    @assert length(a)==14 "a must be a vector of length 14"
    @assert length(b)==14 "b must be a vector of length 14"
end

coordinates_topology(::BosonStarSpacetime) = SphericalTopology()
radius(position, ::BosonStarSpacetime) = position[2]

function metric!(g::AbstractMatrix, position::AbstractVector, spacetime::BosonStarSpacetime)
    r = position[2]
    θ = position[3]

    a = spacetime.a
    b = spacetime.b

    gtt = -1 +
          (1 + r * a[1] + r^2 * a[2] + r^3 * a[3] + r^4 * a[4] + r^5 * a[5] + r^6 * a[6]) /
          (a[7] + r * a[8] + r^2 * a[9] + r^3 * a[10] + r^4 * a[11] + r^5 * a[12] +
           r^6 * a[13] + r^7 * a[14])
    grr = 1 -
          (r * b[1] + r^2 * b[2] + r^3 * b[3] + r^4 * b[4] + r^5 * b[5] + r^6 * b[6]) /
          (b[7] + r * b[8] + r^2 * b[9] + r^3 * b[10] + r^4 * b[11] + r^5 * b[12] +
           r^6 * b[13] + r^7 * b[14])

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

function metric_inverse!(g, position, spacetime::BosonStarSpacetime, gaux, cache)
    r = position[2]
    θ = position[3]

    a = spacetime.a
    b = spacetime.b

    gtt = -1 +
          (1 + r * a[1] + r^2 * a[2] + r^3 * a[3] + r^4 * a[4] + r^5 * a[5] + r^6 * a[6]) /
          (a[7] + r * a[8] + r^2 * a[9] + r^3 * a[10] + r^4 * a[11] + r^5 * a[12] +
           r^6 * a[13] + r^7 * a[14])
    grr = 1 -
          (r * b[1] + r^2 * b[2] + r^3 * b[3] + r^4 * b[4] + r^5 * b[5] + r^6 * b[6]) /
          (b[7] + r * b[8] + r^2 * b[9] + r^3 * b[10] + r^4 * b[11] + r^5 * b[12] +
           r^6 * b[13] + r^7 * b[14])

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

allocate_christoffel_cache(::BosonStarSpacetime) = nothing

function christoffel!(Γ::AbstractArray, position::AbstractVector, spacetime::BosonStarSpacetime)
    #Spacetime coordinates
    r = position[2]
    θ = position[3]

    a = spacetime.a
    b = spacetime.b

    numa = 1 + r * a[1] + r^2 * a[2] + r^3 * a[3] + r^4 * a[4] + r^5 * a[5] + r^6 * a[6]
    dena = a[7] + r * a[8] + r^2 * a[9] + r^3 * a[10] + r^4 * a[11] + r^5 * a[12] +
           r^6 * a[13] + r^7 * a[14]

    numb = r * b[1] + r^2 * b[2] + r^3 * b[3] + r^4 * b[4] + r^5 * b[5] + r^6 * b[6]
    denb = b[7] + r * b[8] + r^2 * b[9] + r^3 * b[10] + r^4 * b[11] + r^5 * b[12] +
           r^6 * b[13] + r^7 * b[14]

    gtt = -1 + numa / dena
    grr = 1 - numb / denb

    ∂r_numa = a[1] + 2r * a[2] + 3r^2 * a[3] + 4r^3 * a[4] + 5r^4 * a[5] + 6r^5 * a[6]
    ∂r_dena = a[8] + 2r * a[9] + 3r^2 * a[10] + 4r^3 * a[11] + 5r^4 * a[12] + 6r^5 * a[13] +
              7r^6 * a[14]

    ∂r_numb = b[1] + 2r * b[2] + 3r^2 * b[3] + 4r^3 * b[4] + 5r^4 * b[5] + 6r^5 * b[6]
    ∂r_denb = b[8] + 2r * b[9] + 3r^2 * b[10] + 4r^3 * b[11] + 5r^4 * b[12] + 6r^5 * b[13] +
              7r^6 * b[14]

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
    spacetime::BosonStarSpacetime,
    rotation_sense)
    #Spacetime coordinates
    r = position[2]

    a = spacetime.a

    numa = 1 + r * a[1] + r^2 * a[2] + r^3 * a[3] + r^4 * a[4] + r^5 * a[5] + r^6 * a[6]
    dena = a[7] + r * a[8] + r^2 * a[9] + r^3 * a[10] + r^4 * a[11] + r^5 * a[12] +
           r^6 * a[13] + r^7 * a[14]

    ∂r_numa = a[1] + 2r * a[2] + 3r^2 * a[3] + 4r^3 * a[4] + 5r^4 * a[5] + 6r^5 * a[6]
    ∂r_dena = a[8] + 2r * a[9] + 3r^2 * a[10] + 4r^3 * a[11] + 5r^4 * a[12] + 6r^5 * a[13] +
              7r^6 * a[14]

    ∂r_gtt = ∂r_numa / dena - numa * ∂r_dena / dena^2
    ∂r_gφφ = 2r

    Ω = sqrt(-∂r_gtt / ∂r_gφφ)

    return sign(rotation_sense) * Ω
end

function BosonStarSpacetime(name::Symbol)
    if name==:LBS1
        a = [-0.12835651132329287, 0.013682816437869968, -0.0000481174714805307,
            0.00005645401289538938, -4.950385183447073e-6, 1.2257668397193955e-6,
            2.234495931121217, -0.286705001384064, 0.0535085065184157,
            -0.002796436948776261, 0.00043953642633610485, 1.5216185242475638e-5,
            -2.2791584430232886e-6, 6.116664949485882e-7]

        b = [0.002169200852262319, -0.19546743355115126, 0.06287130575728342,
            -0.010651955590611452, 0.0009203500154294379, -3.751481658539803e-5,
            13.99470385461245, -3.910737098471579, 0.3926444282051113,
            0.1007509282796355, -0.03817380739469051, 0.00609572731408874,
            -0.0004946554674560391, 1.8733095859225056e-5]
    elseif name==:LBS2
        a = [
            -0.1435113714814284,
            0.017087213017549247,
            0.00006572667868294113,
            0.00006653559130938554,
            -4.06334402727389e-6,
            2.3329241292681612e-6,
            1.9664515659647512,
            -0.28214639097999417,
            0.05920959472144834,
            -0.003337868600178898,
            0.0006589237037906891,
            1.4245462881104672e-5,
            -1.7250207397416203e-6,
            1.1644482387750522e-6,
        ]

        b = [
            0.16195314443269665,
            -25.820415493731403,
            9.333360179608434,
            -1.962592139922535,
            0.21138049431066872,
            -0.01104297071863971,
            1205.626774641147,
            -389.135691872525,
            56.177826933618626,
            10.504677859429323,
            -5.645148389357582,
            1.1476144411789206,
            -0.11570280600991613,
            0.005512416121772507,
        ]

    elseif name==:LBS3

        a = [-0.18333074964802729, 0.02824076525546771, 0.0000630617866798556,
            0.00025089294300200305, -0.000037473999198027776, 0.000012657133473364083,
            1.561206056558541, -0.2860839971711376, 0.07471646699977556,
            -0.004880895957362787, 0.0012816989648470412, 0.00007862954431656809,
            -0.000017823667123517417, 6.321472556491356e-6]

        b = [0.015465297773050266, -2.456296533773868, 1.0463319819861836,
            -0.27891286253840036, 0.03845038475111201, -0.0026622868988117583,
            48.033264323902955, -17.61227848527583, 2.5067857309445802,
            1.2858098202918993, -0.6919234939481852, 0.1714620732009562,
            -0.021705027093548865, 0.0013292192166510165]

    elseif name==:SBS1

        a = [-0.4816713170245124, 0.09522533212597596, -0.011851078458246697,
            0.0018333581130592178, -0.00030102180046519137, 0.000022669117313507083,
            1.7801329298094903, -0.8540698549417561, 0.1792631032506714,
            -0.020894876676239463, 0.00036248443480244735, 0.0005956260834981064,
            -0.00014219669566471017, 0.000011251827048852713]

        b = [3.464440462907006, -77.43547621253862, 52.18873402574934,
            -14.234709268092356, 1.83070437653844, -0.09248708434753089,
            10563.532733322889, -7371.646735710398, 2095.241360836071,
            -257.51183658547023, -7.878431286777774, 7.214853143885797,
            -0.9596134119818048, 0.045718880775679406]

    elseif name==:SBS2

        a = [
            -0.770122726874973,
            0.2429199224348993,
            -0.049258805846667024,
            0.012246484883771286,
            -0.0030531459539308056,
            0.000345616329428982,
            1.2429197582137703,
            -0.9524553681667818,
            0.30714796306515635,
            -0.051135447355958964,
            -0.0006628533482589004,
            0.00430763937639871,
            -0.0014591131511188882,
            0.000171898015431344,
        ]

        b = [
            -7.778311830017875,
            -73.755465597299,
            85.06142818335147,
            -36.062426601741485,
            6.890374421385985,
            -0.5029387517945284,
            2532.806075630009,
            -2991.777460224809,
            1389.163932366156,
            -258.8573147429205,
            -21.934080115129785,
            20.195069470771628,
            -3.748875278619517,
            0.2483481698219909,
        ]

    elseif name==:SBS3

        a = [
            -0.9006548769780998,
            0.33564557693094266,
            -0.08076401411508972,
            0.02443773754042898,
            -0.007634296546364588,
            0.0010971769721096316,
            1.0668212434807347,
            -0.9566244971860363,
            0.3592571346542265,
            -0.06683580469209711,
            -0.0026951710909715577,
            0.008943035421722353,
            -0.0036793009818982428,
            0.0005465290841353206,
        ]
        b = [
            -1.2425056962320251,
            -6.960345160469196,
            10.693372623349484,
            -5.715235155623435,
            1.3547338843187757,
            -0.12160721643561871,
            73.27385618453417,
            -104.16873907198386,
            52.37625328612954,
            -3.1121062971099924,
            -8.323536963898565,
            3.995574425107983,
            -0.7926378556435706,
            0.060820067311090256,
        ]
    else
        error("Unknown name $name")
    end
    return BosonStarSpacetime(a, b)
end
