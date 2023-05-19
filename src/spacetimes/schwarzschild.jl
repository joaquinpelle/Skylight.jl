#Kerr-Schild coordinates

@with_kw struct SchwarzschildSpacetimeKerrSchildCoordinates <: AbstractSpacetime 

    M::Float64
    @assert M >= 0.0

    #Cache
    l::Vector{Float64}=zeros(4)

end

coordinate_system_class(::SchwarzschildSpacetimeKerrSchildCoordinates) = CartesianClass()

function set_metric!(g, position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates)

    """ 
    g: container for the metric 
    q: spacetime position
    """

    M = spacetime.M
    
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    H = 2M/r
    
    l = spacetime.l

    l[1] = 1.0
    l[2] = x/r
    l[3] = y/r
    l[4] = z/r
    
    g[1,1]=-1. + H*l[1]*l[1]
    g[1,2]= 0. + H*l[1]*l[2]
    g[1,3]= 0. + H*l[1]*l[3]
    g[1,4]= 0. + H*l[1]*l[4]
    g[2,1]= g[1,2]
    g[2,2]= 1. + H*l[2]*l[2]
    g[2,3]= 0. + H*l[2]*l[3]
    g[2,4]= 0. + H*l[2]*l[4]
    g[3,1]= g[1,3]
    g[3,2]= g[2,3]
    g[3,3]= 1. + H*l[3]*l[3]
    g[3,4]= 0. + H*l[3]*l[4]
    g[4,1]= g[1,4]
    g[4,2]= g[2,4]
    g[4,3]= g[3,4]
    g[4,4]= 1. + H*l[4]*l[4]
    
    return nothing
    
end

function set_metric_inverse!(g, position, spacetime::SchwarzschildSpacetimeKerrSchildCoordinates)

    """ 
    g: container for the metric 
    q: spacetime position
    """

    M = spacetime.M
    
    t, x, y, z = position
    r = sqrt(x^2 + y^2 + z^2)
    H = 2M/r
    
    l = spacetime.l
    
    l[1] =-1.
    l[2] = x/r
    l[3] = y/r
    l[4] = z/r
    
    g[1,1]=-1. - H
    g[1,2]= 0. - H*l[1]*l[2]
    g[1,3]= 0. - H*l[1]*l[3]
    g[1,4]= 0. - H*l[1]*l[4]
    g[2,1]= g[1,2]
    g[2,2]= 1. - H*l[2]*l[2]
    g[2,3]= 0. - H*l[2]*l[3]
    g[2,4]= 0. - H*l[2]*l[4]
    g[3,1]= g[1,3]
    g[3,2]= g[2,3]
    g[3,3]= 1. - H*l[3]*l[3]
    g[3,4]= 0. - H*l[3]*l[4]
    g[4,1]= g[1,4]
    g[4,2]= g[2,4]
    g[4,3]= g[3,4]
    g[4,4]= 1. - H*l[4]*l[4]
    
    return nothing
    
end

event_horizon_radius(spacetime::SchwarzschildSpacetimeKerrSchildCoordinates) = 2*spacetime.M

#Spherical coordinates

@with_kw struct SchwarzschildSpacetimeSphericalCoordinates <: AbstractSpacetime 

    M::Float64
    @assert M >= 0.0

end

coordinate_system_class(::SchwarzschildSpacetimeSphericalCoordinates) = SphericalClass()

function set_metric!(g, q, spacetime::SchwarzschildSpacetimeSphericalCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """

    t, r, θ, φ = q

    M = spacetime.M

    fill!(g,0.0)
    g[1,1] = -(1-2M/r)
    g[2,2] =  1/(1-2M/r)
    g[3,3] =  r^2
    g[4,4] =  r^2*sin(θ)^2
    
    return nothing

end

function set_metric_inverse!(g, q, spacetime::SchwarzschildSpacetimeSphericalCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """

    t, r, θ, φ = q

    M = spacetime.M

    fill!(g,0.0)
    g[1,1] = -1/(1-2M/r)
    g[2,2] =  1-2M/r
    g[3,3] =  1/r^2
    g[4,4] =  1/(r^2*sin(θ)^2)
    
    return nothing

end

event_horizon_radius(spacetime::SchwarzschildSpacetimeSphericalCoordinates) = 2*spacetime.M
