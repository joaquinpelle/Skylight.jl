# Automatic differentiation

Unless a more specialized method is available, the Christoffel symbols of a given spacetime
are calculated using forward-mode automatic differentiation from the metric coefficients, using
functionality from the ForwardDiff package. 

```@docs
Skylight.AutoDiffChristoffelCache
```