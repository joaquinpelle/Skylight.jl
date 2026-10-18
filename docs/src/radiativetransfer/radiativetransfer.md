# Radiative transfer
## Physical units in volume transfer

Volume emissivities are per photon energy in CGS,
`j_E` in erg s⁻¹ cm⁻³ sr⁻¹ erg⁻¹, and absorptivities are in cm⁻¹.
`observation_energies` are in erg, while coordinates and the geodesic affine
parameter use the geometrized length unit selected by
`unit_mass_in_solar_masses`. The initialized past-directed momentum has unit
energy in the camera tetrad (or the reference frame used by the image plane).

For `q = u ⋅ k`, `E = E0*q`, and `L = G*M_unit/c²` in cm, the path length measured
by the plasma is `dl = L*q*ds`. The volume state stores optical depth and the
accumulated invariant intensity `J = I_E/E³`. Its derivatives are

```math
\frac{d\tau}{ds}=\frac{L}{E_0}E\alpha_E,\qquad
\frac{dJ}{ds}=\frac{L}{E_0}\frac{j_E}{E^2}e^{-\tau}.
```

The postprocessor restores the observer's energy cubed. Converting a spectrum
per energy to one per frequency requires `I_ν = h*I_E` and `F_ν = h*F_E`.
Do not add another length conversion to CGS coefficients in a model.

Volume integration has no optical-depth cutoff by default (`τmax=Inf`). A finite
`τmax` terminates a ray only after every frequency exceeds it. Choose and check
this approximation against the source function and the desired accuracy;
`exp(-τmax)` alone is not an absolute error bound.
