#TODO add reference

@with_kw struct BosonStarSpacetime <: AbstractRegularCompactObjectSpacetime
    a::Vector{Float64}
    b::Vector{Float64}

    @assert length(a) == 14 "a must be a vector of length 14"
    @assert length(b) == 14 "b must be a vector of length 14"
end

coordinates_topology(::BosonStarSpacetime) = SphericalTopology() 
radius(position, ::BosonStarSpacetime) = position[2]

function metric!(g, position, spacetime::BosonStarSpacetime)
    t, r, θ, φ = position

    a = spacetime.a
    b = spacetime.b

    gtt = -1+(1+r*a[1]+r^2*a[2]+r^3*a[3]+r^4*a[4]+r^5*a[5]+r^6*a[6])/(a[7]+r*a[8]+r^2*a[9]+r^3*a[10]+r^4*a[11]+r^5*a[12]+r^6*a[13]+r^7*a[14])
    grr =  1-(r*b[1]+r^2*b[2]+r^3*b[3]+r^4*b[4]+r^5*b[5]+r^6*b[6])/(b[7]+r*b[8]+r^2*b[9]+r^3*b[10]+r^4*b[11]+r^5*b[12]+r^6*b[13]+r^7*b[14])
    
    g[1,1]= gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= grr 
    g[2,3]= 0.
    g[2,4]= 0.
    g[3,1]= 0.
    g[3,2]= 0.
    g[3,3]= r^2
    g[3,4]= 0.
    g[4,1]= 0.
    g[4,2]= 0.
    g[4,3]= 0.
    g[4,4]= r^2*sin(θ)^2
    
    return nothing
end

function metric_inverse!(g, position, spacetime::BosonStarSpacetime)
    t, r, θ, φ = position

    a = spacetime.a
    b = spacetime.b

    gtt = -1+(1+r*a[1]+r^2*a[2]+r^3*a[3]+r^4*a[4]+r^5*a[5]+r^6*a[6])/(a[7]+r*a[8]+r^2*a[9]+r^3*a[10]+r^4*a[11]+r^5*a[12]+r^6*a[13]+r^7*a[14])
    grr =  1-(r*b[1]+r^2*b[2]+r^3*b[3]+r^4*b[4]+r^5*b[5]+r^6*b[6])/(b[7]+r*b[8]+r^2*b[9]+r^3*b[10]+r^4*b[11]+r^5*b[12]+r^6*b[13]+r^7*b[14])
    
    g[1,1]= 1/gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= 1/grr 
    g[2,3]= 0.
    g[2,4]= 0.
    g[3,1]= 0.
    g[3,2]= 0.
    g[3,3]= 1.0/r^2
    g[3,4]= 0.
    g[4,1]= 0.
    g[4,2]= 0.
    g[4,3]= 0.
    g[4,4]= 1.0/(r^2*sin(θ)^2)
    
    return nothing
end

allocate_christoffel_cache(::BosonStarSpacetime) = nothing

function christoffel!(Γ, position, spacetime::BosonStarSpacetime)
    #Spacetime coordinates
    t, r, θ, φ = position

    a = spacetime.a
    b = spacetime.b

    numa = 1+r*a[1]+r^2*a[2]+r^3*a[3]+r^4*a[4]+r^5*a[5]+r^6*a[6]
    dena = a[7]+r*a[8]+r^2*a[9]+r^3*a[10]+r^4*a[11]+r^5*a[12]+r^6*a[13]+r^7*a[14]
    
    numb = r*b[1]+r^2*b[2]+r^3*b[3]+r^4*b[4]+r^5*b[5]+r^6*b[6]
    denb = b[7]+r*b[8]+r^2*b[9]+r^3*b[10]+r^4*b[11]+r^5*b[12]+r^6*b[13]+r^7*b[14]
 
    gtt = -1+numa/dena
    grr =  1-numb/denb

    ∂r_numa = a[1]+2r*a[2]+3r^2*a[3]+4r^3*a[4]+5r^4*a[5]+6r^5*a[6]
    ∂r_dena = a[8]+2r*a[9]+3r^2*a[10]+4r^3*a[11]+5r^4*a[12]+6r^5*a[13]+7r^6*a[14]
    
    ∂r_numb = b[1]+2r*b[2]+3r^2*b[3]+4r^3*b[4]+5r^4*b[5]+6r^5*b[6]
    ∂r_denb = b[8]+2r*b[9]+3r^2*b[10]+4r^3*b[11]+5r^4*b[12]+6r^5*b[13]+7r^6*b[14]

    ∂r_gtt =  ∂r_numa/dena - numa*∂r_dena/dena^2
    ∂r_grr = -∂r_numb/denb + numb*∂r_denb/denb^2

    dα = ∂r_gtt/(2*gtt)
    dβ = ∂r_grr/(2*grr)
    
    Γ[1,1,2] = dα
    Γ[1,2,1] = Γ[1,1,2]
    
    Γ[2,1,1] = -gtt/grr*dα
    Γ[2,2,2] =  dβ
    Γ[2,3,3] = -r/grr
    Γ[2,4,4] = -r/grr*sin(θ)^2
    
    Γ[3,2,3] =  1.0/r
    Γ[3,4,4] = -sin(θ)*cos(θ)
    Γ[3,3,2] =  Γ[3,2,3]
    
    Γ[4,2,4] = 1.0/r
    Γ[4,3,4] = cos(θ)/sin(θ)
    Γ[4,4,2] = Γ[4,2,4]
    Γ[4,4,3] = Γ[4,3,4]
    
    return nothing
end

function circular_geodesic_angular_speed(position, spacetime::BosonStarSpacetime, rotation_sense)
    #Spacetime coordinates
    t, r, θ, φ = position

    a = spacetime.a

    numa = 1+r*a[1]+r^2*a[2]+r^3*a[3]+r^4*a[4]+r^5*a[5]+r^6*a[6]
    dena = a[7]+r*a[8]+r^2*a[9]+r^3*a[10]+r^4*a[11]+r^5*a[12]+r^6*a[13]+r^7*a[14]
    
    ∂r_numa = a[1]+2r*a[2]+3r^2*a[3]+4r^3*a[4]+5r^4*a[5]+6r^5*a[6]
    ∂r_dena = a[8]+2r*a[9]+3r^2*a[10]+4r^3*a[11]+5r^4*a[12]+6r^5*a[13]+7r^6*a[14]
    
    ∂r_gtt =  ∂r_numa/dena - numa*∂r_dena/dena^2
    ∂r_gφφ = 2r

    Ω = sqrt(-∂r_gtt/∂r_gφφ)

    return sign(rotation_sense)*Ω
end