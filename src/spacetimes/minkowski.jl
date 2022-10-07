export MinkowskiSpacetimeCartesianCoordinates
export MinkowskiSpacetimeSphericalCoordinates

struct MinkowskiSpacetimeCartesianCoordinates <: AnalyticSpacetime end

coordinate_system_kind(spacetime::MinkowskiSpacetimeCartesianCoordinates) = CartesianKind()

function set_metric!(g, q, spacetime::MinkowskiSpacetimeCartesianCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """

    fill!(g,0.0)
    g[1,1] = -1.0
    g[2,2] =  1.0
    g[3,3] =  1.0
    g[4,4] =  1.0

    return nothing
    
end

set_metric_inverse!(g, q, spacetime::MinkowskiSpacetimeCartesianCoordinates) = set_metric!(g, q, spacetime)

#Spherical coordinates

struct MinkowskiSpacetimeSphericalCoordinates <: AnalyticSpacetime end

coordinate_system_kind(spacetime::MinkowskiSpacetimeSphericalCoordinates) = SphericalKind()

function set_metric!(g, q, spacetime::MinkowskiSpacetimeSphericalCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """

    t, r, θ, φ = q

    fill!(g,0.0)
    g[1,1] = -1.0
    g[2,2] =  1.0
    g[3,3] =  r^2
    g[4,4] =  r^2*sin(θ)^2
    
    return nothing

end

function set_metric_inverse!(g, q, spacetime::MinkowskiSpacetimeSphericalCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """

    t, r, θ, φ = q

    fill!(g,0.0)
    g[1,1] = -1.0
    g[2,2] =  1.0
    g[3,3] =  1/r^2
    g[4,4] =  1/(r^2*sin(θ)^2)
    
    return nothing

end