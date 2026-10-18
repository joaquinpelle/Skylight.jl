# Four ion-torus configurations

`run/iontorus.jl` runs the following controlled comparison, using `M=1` in code units:

| Case | a/M | R0 M² | cusp/M | center/M | outer/M | central electron temperature (K) |
|---|---:|---:|---:|---:|---:|---:|
| C1 | 0.10 | 0 | 4.01007 | 8.92492 | 34.91644 | 6.83397e8 |
| C2 | 0.10 | -0.001 | 4.01895 | 8.27045 | 18.32472 | 6.10734e8 |
| C3 | 0.99 | 0 | 1.23674 | 1.87498 | 5.81916 | 3.29988e9 |
| C4 | 0.99 | -0.001 | 1.23548 | 1.87048 | 5.71326 | 3.29841e9 |

All cases use prograde rotation, `lambda=0.7`, central mass density `1e-17 g/cm³`,
polytropic index `n=1.5`, magnetic pressure fraction `beta=0.5`, central
temperature ratio `xi=0.01`, hydrogen/helium mass fractions 0.75/0.25, and a
mass unit of `4.14e6` solar masses. Emission includes thermal synchrotron and
bremsstrahlung, with thermal self-absorption. There is no Compton scattering.

The static pinhole camera is at `r=50 M`, inclination 60 degrees, with full
angular apertures 120 by 90 degrees. Frequencies measured by that observer are
86, 230 and 345 GHz. This is a finite-distance comparison: the flux is local to
the camera, especially for the asymptotically AdS cases. It is not an Earth flux
or an angular size in microarcseconds.

## Running

From the repository root, install the project dependencies once:

```sh
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

Validate parameters, run a preview, or select a case:

```sh
julia --project=. run/iontorus.jl --dry-run
julia --project=. --threads=auto run/iontorus.jl --case all --pixels 256 --output results/preview
julia --project=. --threads=auto run/iontorus.jl --case C4 --pixels 1024 --output results/production
```

The horizontal pixel count must be divisible by four; the vertical count is
three quarters of it. Use `--case all --pixels 1024` for all four production
images. Production runs can be expensive. `--help` lists all options. Including
the file from Julia defines the `IonTorusExample` module without starting a run.
CairoMakie is not required.

Each case writes an HDF5 file, a CSV containing the three frequency samples,
and a grayscale SVG preview at 230 GHz. Existing results are protected: choose
a new output directory for another run. HDF5 stores `I_E`, `I_nu`, optical depth,
local `F_nu`, the initial/final ray states, angular axes, camera tetrad and model
and solver parameters. Image arrays are indexed as `(frequency, alpha, beta)`.
Units are stored as file attributes. `I_nu = h*I_E`; the accumulated invariant
intensity in the ODE state is not itself an observed intensity. The SVG uses a
six-decade logarithmic stretch without modifying the scientific arrays.

## Numerical checks

The default step limit is `0.25 M`, with `reltol=1e-8`, `abstol=1e-10`, and
termination `0.001 M` above the event horizon. The script verifies that this
cutoff is below the torus cusp. The previous `0.2 M` cutoff removes part of the
high-spin tori. Optical depth has no termination cutoff (`tau_max=Inf`). Failed
ODE solves raise an error instead of being silently accepted.

The intensity components use a separate absolute tolerance, specified in
observable units by `--intensity-abstol 1e-18` (CGS specific intensity per Hz).
The script converts this to the corresponding tolerance on `I_E/E³`. A single
scalar tolerance for all components can stall initially dark rays at the torus
boundary. Finite radio invariant intensities can exceed `1e50`; the nonvacuum
solver checks finiteness without the dimensionful magnitude cutoff used by
some older DiffEq releases.

Before interpreting production differences, repeat with twice the angular
resolution, `--dtmax 0.125`, `--horizon-buffer 0.0005`, and tighter tolerances,
and compare the saved fluxes and resolved features. The regression suite checks
analytic transfer, equilibrium, geometry and small end-to-end images; it does
not establish convergence of a production image.

## Physical conventions and limits

`IonTorus(spacetime)` now enforces the integrated hydrostatic equation with
`e ≈ rho*c²` and a total-pressure polytrope. The central temperature follows
from the potential depth. To reproduce an independently prescribed temperature,
request `thermodynamics=:prescribed, Tec=...` explicitly; that profile does not
generally satisfy hydrostatic balance. The raw `IonTorus(;...)` constructor has
no geometry and remains a parameter container in prescribed mode.

For `R0<=0`, the angular-momentum reference uses the **inner** circular-orbit
root of the rescaled Killing energy `E=-Xi*u_t=1`. In AdS this is a convention,
not a marginal escape orbit. Positive-curvature cusp-filled tori are rejected:
they need the outer potential barrier as well, which is not implemented here.
Degenerate horizons are outside the nonzero-curvature root selection used here.

The legacy synchrotron fit neglects emission below `5e8 K`. It interpolates the
tabulated coefficients up to `3e10 K`, joins them linearly to the relativistic
limit at `3.2e10 K`, and uses that limit above it. The join is an explicit
numerical prescription. In particular, the cold cutoff can affect the outer
parts of these tori and must be assessed before quantitative interpretation.
Bremsstrahlung includes the charge-squared weighted ion density for fully
ionized H/He. Thermal absorptivities obey Kirchhoff's law.

The equilibrium construction follows [Straub et al. (2012)](https://www.aanda.org/articles/aa/abs/2012/07/aa19209-12/aa19209-12.html).
