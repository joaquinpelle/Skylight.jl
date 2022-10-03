export MinkowskiSpacetimeParameters 
export MinkowskiSpacetimeCartesianCoordinates, MinkowskiSpacetimeSphericalCoordinates

struct MinkowskiSpacetimeParameters <: SpacetimeParameters end

#Cartesian coordinates

@with_kw struct MinkowskiSpacetimeCartesianCoordinates{T<:Function} <: AnalyticSpacetime

    parameters::MinkowskiSpacetimeParameters = MinkowskiSpacetimeParameters() 
    coordinate_system_kind::CartesianKind = CartesianKind()
    metric!::T = minkowski_metric_cartesian_coordinates!

end

function minkowski_metric_cartesian_coordinates!(g, q, par::MinkowskiSpacetimeParameters)
    
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

@with_kw struct MinkowskiSpacetimeSphericalCoordinates{T<:Function} <: AnalyticSpacetime 

    parameters::MinkowskiSpacetimeParameters = MinkowskiSpacetimeParameters()
    coordinate_system_kind::SphericalKind = SphericalKind()
    metric!::T = minkowski_metric_spherical_coordinates!

end

function minkowski_metric_spherical_coordinates!(g, q, par::MinkowskiSpacetimeParameters)
    
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