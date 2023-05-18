@with_kw struct StarAcrossWormholeCallbackParameters <: AbstractCallbackParameters
    
    rmax::Float64
    l_center::Float64
    star_radius::Float64

end

function get_cb_params(model::StarAcrossWormhole, configurations) 

    l_center = model.l_center
    star_radius = model.star_radius
    rmax = get_rmax(configurations)

    return StarAcrossWormholeCallbackParameters(l_center=l_center, star_radius=star_radius, rmax=rmax)

end

get_callback(::StarAcrossWormhole, ) = star_across_wormhole_callback()

star_across_wormhole_callback() = VectorContinuousCallback(star_across_wormhole_condition, star_across_wormhole_affect!, 2)

function star_across_wormhole_condition(out, u, t, integrator)
    
    l_center = integrator.p.cb_params.l_center
    star_radius = integrator.p.cb_params.star_radius
    rmax = integrator.p.cb_params.rmax

    b0 = integrator.p.spacetime.b0
    Q = integrator.p.spacetime.Q

    t,l,θ,φ = u[1:4]

    r = sqrt(l^2+b0^2-Q^2)

    x,y,z = cartesian_from_spherical([r,θ,φ])
    
    r_center = sqrt(l_center^2+b0^2-Q^2)
    
    out[1] = (x-r_center)^2 + y^2 + z^2 - star_radius^2  
    out[2] = (rmax-r)*(rmax+r)

end

function star_across_wormhole_affect!(integrator,idx)
 
    if idx==1 && integrator.u[2]<0.0  
        terminate!(integrator)
    elseif idx==2 
        terminate!(integrator) 
    end

end
