# Skylight.jl

Skylight is a Julia package for general-relativistic ray-tracing and radiative transfer in curved spacetimes. It supports arbitrary spacetimes, not relying on any symmtries nor asymptotic flatness. It is designed to be fast, accurate and easily extensible to user-defined spacetimes and radiative models. It uses [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) to compute the Christoffel symbols from the metric, and supports multithreading in its compute-intensive routines. At its core, Skylight simultaneously solves the geodesic equations and the covariant transport equations along the geodesics, i.e. 

```math
\frac{\mathop{d^2 x^{\alpha}}}{\mathop{d\lambda^2}}+\Gamma^\alpha_{\mu \nu} \frac{\mathop{d x^{\mu}}}{\mathop{d\lambda}} \frac{\mathop{d x^{\nu}}}{\mathop{d\lambda}}=0\,, 
```

where $\Gamma^\alpha_{\mu \nu}$ are the Christoffel symbols of the spacetime, $x^\alpha$ is the position and $\lambda$ is an affine parameter along the geodesic, together with

```math
\frac{\mathop{d}}{\mathop{d\lambda}} \left( \frac{I_\nu}{\nu^3}\right) = \frac{j_\nu}{\nu^2} - \nu \alpha_\nu \left( \frac{I_\nu}{\nu^3}\right) \,,
```

where $\nu$ is the frequency, $I_\nu$ is the intensity of the radiation field, and $j_\nu$ and $\alpha_\nu$ are the emissivity and absorptivity coefficients of the medium, respectively. Skylight also has special functions for surface emission models with transport in vacuum, where the transport can be much simplified to just connecting the intensity of the radiation field between the emission and observation points.     

See [Getting started](@ref) for a quick start guide.

This documentation is under construction.