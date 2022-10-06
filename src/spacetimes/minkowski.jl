export MinkowskiSpacetimeParameters 
export MinkowskiSpacetimeCartesianCoordinates, MinkowskiSpacetimeSphericalCoordinates

struct MinkowskiSpacetimeParameters <: SpacetimeParameters end

#Cartesian coordinates

@with_kw struct MinkowskiSpacetimeCartesianCoordinates <: AnalyticSpacetime

    parameters::MinkowskiSpacetimeParameters = MinkowskiSpacetimeParameters() 
    coordinate_system_kind::CartesianKind = CartesianKind()

end

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

@with_kw struct MinkowskiSpacetimeSphericalCoordinates <: AnalyticSpacetime 

    parameters::MinkowskiSpacetimeParameters = MinkowskiSpacetimeParameters()
    coordinate_system_kind::SphericalKind = SphericalKind()

end

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