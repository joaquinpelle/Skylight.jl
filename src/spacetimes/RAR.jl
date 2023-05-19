@with_kw struct RARInterpolator{T}
    
    gtt :: T
    grr :: T
    dν :: T
    M :: T
    dM :: T

end

function RARInterpolator(data_dir)
        
    #Data files must be in geometrized units
    gtt = readdlm(data_dir*"/gtt.txt", ' ', Float64, '\n')
    grr = readdlm(data_dir*"/grr.txt", ' ', Float64, '\n')
    dν = readdlm(data_dir*"/dnu.txt", ' ', Float64, '\n')
    M = readdlm(data_dir*"/M.txt", ' ', Float64, '\n')
    dM = readdlm(data_dir*"/dM.txt", ' ', Float64, '\n')
    
    return RARInterpolator(gtt=CubicSpline(-gtt[:,2], gtt[:,1]), #with the minus sign for the (- + + +) convention
                           grr = CubicSpline(grr[:,2], grr[:,1]),
                           dν = CubicSpline(dν[:,2], dν[:,1]),
                           M = CubicSpline(M[:,2], M[:,1]),
                           dM = CubicSpline(dM[:,2], dM[:,1]))

end

@with_kw struct RARSpacetime{T} <: AbstractSpacetime
    
    data_dir::String  
    interp::RARInterpolator{T} = RARInterpolator(data_dir)

end

coordinate_system_class(::RARSpacetime) = SphericalClass()

function set_metric!(g,point, spacetime::RARSpacetime)
    
    """"Ruffini-Arguelles-Rueda metric for dark-matter halo"""
    
    t, r, θ, φ = point
    
    interp = spacetime.interp

    gtt = interp.gtt(r)
    grr = interp.grr(r)
    
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

function set_metric_inverse!(g, point, spacetime::RARSpacetime)
    
    t, r, θ, φ = point

    interp = spacetime.interp
    
    gtt = interp.gtt(r)
    grr = interp.grr(r)
    
    g[1,1]= 1.0/gtt
    g[1,2]= 0. 
    g[1,3]= 0. 
    g[1,4]= 0. 
    g[2,1]= 0.
    g[2,2]= 1.0/grr 
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

function set_christoffel!(Γ, position, spacetime::RARSpacetime)

    t, r, θ, φ = position
    
    interp = spacetime.interp

    gtt = interp.gtt(r)
    grr = interp.grr(r)
    dν = interp.dν(r)
    M = interp.M(r)
    dM = interp.dM(r)
    
    dλ = 2*grr/r*(dM-M/r)
    
    Γ[1,1,2] = 0.5*dν
    Γ[1,2,1] = Γ[1,1,2]
    
    Γ[2,1,1] = -gtt/grr*0.5*dν
    Γ[2,2,2] = 0.5*dλ
    Γ[2,3,3] = -r/grr
    Γ[2,4,4] = -r/grr*sin(θ)^2
    
    Γ[3,2,3] = 1.0/r
    Γ[3,4,4] = -sin(θ)*cos(θ)
    Γ[3,3,2] = Γ[3,2,3]
    
    Γ[4,2,4] = 1.0/r
    Γ[4,3,4] = cos(θ)/sin(θ)
    Γ[4,4,2] = Γ[4,2,4]
    Γ[4,4,3] = Γ[4,3,4]
    
    return nothing

end


