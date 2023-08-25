In the observer-to-emitter scheme the initial momenta are past-directed and pointing inwards. This is valid because of the following reason. The Lioville vector field on the tangent bundle is the vector field which generates the geodesic flow. The component of this field along the fibers is invariant to sign inversion on the fiber. This means
that the geodesic passing through the point $(x^μ,k^μ)$ and the one passing through $(x^μ,-k^μ)$ project to the same curve on the base spacetime. We do it this way so the geodesic integrator kernel is common to both transport schemes. Otherwise, although the physical equations would be the same, the numerical integrators would be different because we would have to take backward steps in the observer-to-emitter case. Once we have the solutions, we invert the sign of $k^μ$ when we use it elsewhere. 
is_position_at_source should't be used in condition/affect. Its purpose is to classify the endstate
of a geodesic already integrated, not to condition the geodesic integration.

Currently Christoffels assume the input array is filled with zeros

Callback params needs to have an rmax. Fix or make explicit.

Coordinates are assumed to be have either cartesian or spherical topology, and to be ordered like
t, r, θ, φ or t, x, y, z. The first coordinate is assumed to be timelike, and the other three spatial.

AbstractAccretionDisk assumes inner_radius and outer_radius as fields.

The procedure is modularized sufficiently so that custom initial data can be provided as well to the integrator (as a (8, nrays) array)

The integrator can integrate geodesics of any kind (timelike and spacelike too)

You need invariant emissivity to work for vectors of energies
Setting absorptovity function to nothing is equivalent to setting absorptivity to zero

The intensity integrated in non-vacuum is the invariant intensity

The non-vacuum transfer equations only work towards the past because of rest frame energy sign assumption

Observation energies in non-vacuum need to be their CGS values