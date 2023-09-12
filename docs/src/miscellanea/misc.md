# Miscellanea

In the observer-to-emitter scheme the initial momenta are past-directed and pointing inwards. To see why this is valid cosider the following: the Lioville vector field on the tangent bundle of spacetime is the vector field which generates the geodesic flow. The component of this field along the fibers is invariant to sign inversion on the fiber. This means that the geodesic passing through the point $(x^μ,k^μ)$ and the one passing through $(x^μ,-k^μ)$ project to the same curve on the base spacetime. We do it this way so the geodesic integrator kernel is common to both transport schemes. Otherwise, although the physical equations would be the same, the numerical integrators would have to be different because we would have to take backward steps in the observer-to-emitter case. Once we have the solutions, we can invert the sign of $k^μ$ if necessary when we use it elsewhere. 

is_position_at_source should't be used in condition/affect. Its purpose is to classify the endstate
of a geodesic already stopped, not to condition the geodesic integration.

Callback params needs to have an rmax.

Coordinates are assumed to be have either cartesian or spherical topology, and to be ordered like
$(t, r, \theta, \varphi)$ or $(t, x, y, z)$. The first coordinate is assumed to be temporal, and the other three spatial. This is checked for in the initialization via the metric signature.

`AbstractAccretionDisk` assumes inner_radius and outer_radius as fields.

The integrator can integrate geodesics of any kind (timelike and spacelike too)

Custom initial data can be provided as well to the integrator as an array of size $(8, N)$ where N is the number of rays.

The intensity integrated in non-vacuum is the invariant intensity

The emissivity must take vectors of energies as inputs

Setting the absorptovity function to return `nothing` is equivalent to setting the absorptivity to zero

The non-vacuum transfer equations only work towards the past because of rest frame energy sign assumption. This will be generalized.

Observation energies in non-vacuum are assumed to be in CGS.