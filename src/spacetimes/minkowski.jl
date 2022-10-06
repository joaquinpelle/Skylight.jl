export MinkowskiSpacetimeCartesianCoordinates
export MinkowskiSpacetimeSphericalCoordinates

struct MinkowskiSpacetimeCartesianCoordinates <: AnalyticSpacetime end

coordinate_system_kind(spacetime::MinkowskiSpacetimeCartesianCoordinates) = CartesianKind()

function set_metric!(g, q, spacetime::MinkowskiSpacetimeCartesianCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """
    @. g = [-1.0 0.0 0.0 0.0;
          0.0 1.0 0.0 0.0;
          0.0 0.0 1.0 0.0;
          0.0 0.0 0.0 1.0]

    return g
    
end

#Spherical coordinates

struct MinkowskiSpacetimeSphericalCoordinates <: AnalyticSpacetime end

coordinate_system_kind(spacetime::MinkowskiSpacetimeSphericalCoordinates) = SphericalKind()

function set_metric!(g, q, spacetime::MinkowskiSpacetimeSphericalCoordinates)
    
    """ 
    g: container for the metric 
    q: spacetime position
    """

    t, r, θ, φ = q

    @. g = [-1.0 0.0 0.0 0.0;
            0.0 1.0 0.0 0.0;
            0.0 0.0 r^2 0.0;
            0.0 0.0 0.0 r^2*sin(θ)^2]
    
    return g

end