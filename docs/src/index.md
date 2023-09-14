# Skylight.jl

(The documentation is under construction...)

Skylight is a Julia package for general-relativistic ray-tracing and radiative transfer in curved spacetimes. It works with any spacetime geometry, without the constraints of specific symmetries or the assumption of asymptotic flatness. It is designed with the following goals in mind:
- Fast computational speed
- High accuracy
- Easy extensibility to user-defined spacetimes and radiative models

It uses [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) from [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) to compute the Christoffel symbols from the spacetime metric, and has multithreading parallelism in its compute-intensive tasks. At its core, Skylight simultaneously solves the geodesic equations and the covariant transport equations along the geodesics, i.e.

```math
\frac{\mathop{d^2 x^{\alpha}}}{\mathop{d\lambda^2}}+\Gamma^\alpha_{\mu \nu} \frac{\mathop{d x^{\mu}}}{\mathop{d\lambda}} \frac{\mathop{d x^{\nu}}}{\mathop{d\lambda}}=0\,, 
```

where $\Gamma^\alpha_{\mu \nu}$ are the Christoffel symbols of the spacetime, $x^\alpha$ is the position and $\lambda$ is an affine parameter along the geodesic, together with

```math
\frac{\mathop{d}}{\mathop{d\lambda}} \left( \frac{I_\nu}{\nu^3}\right) = \frac{j_\nu}{\nu^2} - \nu \alpha_\nu \left( \frac{I_\nu}{\nu^3}\right) \,,
```

where $\nu$ is the frequency, $I_\nu$ is the intensity of the radiation field, and $j_\nu$ and $\alpha_\nu$ are the emissivity and absorptivity coefficients of the medium, respectively. Skylight also has special functions for surface emission models with transport in vacuum, where the transport can be much simplified to just connecting the intensity of the radiation field between the emission and observation points. The inegration of the equations integration is performed with [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). 

See [Getting started](@ref) for a quick start guide.

### Features

#### Supported spacetimes

* Minkowski spacetime
* Schwarzschild spacetime
* Kerr spacetime
* Johannsen spacetime
* f(R)-Kerr spacetime
* Ruffini-Argüelles-Rueda spacetime for fermionic dark matter
* Boson star spacetime with quartic self-interaction and solitonic potentials  
* Extensibility to user-defined spacetimes

#### Radiative models

* Shakura-Sunyaev accretion disks
* Geometrically-thin optically-thick accretion disks with user-provided tabulated temperatures
* Ion torus with synchrotron and bremsstrahlung emission 
* Line emission from accretion disks with user-provided emissivity profiles
* Lamppost corona emission and accretion disk illumination profiles
* Extensibility to user-defined radiative models

#### Geometric and dynamical tools

* Geometric quantities like spacetime metrics, inverse metrics, volume elements, Christoffel symbols, etc.
* Four-vector scalar products, index raising/lowering, orthogonal projection, normalization, etc.
* Constants of motion for spacetimes with symmetries
* Characteristic radii in certain spacetimes, like event horizons, ISCOs, etc. 
* Spacetime geodesics integration

#### Radiative transfer mechanisms

* Radiative transfer in vacuum and in emissive/absorptive media

#### Observable quantities

* Bolometric and specific intensities in generic frames
* Fluxes through arbitrarily oriented surface elements
* Spectra

#### Utilities

* Data loading/saving from/to HDF5 files
* Units and dimensions management

## Installation

#### Requirements
* Julia version at least 1.6

The package is not yet available in the Julia registries. To install it, follow these steps:

1. Clone the repository: `git clone https://github.com/joaquinpelle/Skylight.jl.git`
2. Open the Julia REPL and enter package mode by typing `]`.
3. Add Skylight to your Pkg environment: `] dev \path\to\the\repository`
4. Import Skylight: `using Skylight`

### Folder contents

    .
    ├── docs                    # Documentation files
    ├── run                     # Example scripts and notebooks to run the code
    ├── src                     # Source files
    ├── test                    # Test files 
    │   ├── benchmarks          # Load and stress tests
    │   ├── integration         # End-to-end, integration tests
    │   └── unit                # Unit tests

### To run the unit tests

1. In your terminal, go to the package directory 
2. Open the Julia REPL
3. Go to the package mode by typing `]`
4. Activate the package environment by running the command `activate .`
5. In the package mode, run the command `test`

### To cite this work

If you use this software in your work, we kindly request that you cite [the following paper](https://academic.oup.com/mnras/article-abstract/515/1/1316/6631564)

```bibtex
@article{pelle2022skylight,
  title={Skylight: a new code for general-relativistic ray-tracing and radiative transfer in arbitrary space--times},
  author={Pelle, Joaquin and Reula, Oscar and Carrasco, Federico and Bederian, Carlos},
  journal={Monthly Notices of the Royal Astronomical Society},
  volume={515},
  number={1},
  pages={1316--1327},
  year={2022},
  publisher={Oxford University Press}
}
```