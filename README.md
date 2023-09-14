# Skylight.jl
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joaquinpelle.github.io/Skylight.jl/dev)
[![Build Status](https://github.com/joaquinpelle/Skylight.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/joaquinpelle/Skylight.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/joaquinpelle/Skylight.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/joaquinpelle/Skylight.jl)
[![Coverage](https://coveralls.io/repos/github/joaquinpelle/Skylight.jl/badge.svg?branch=main)](https://coveralls.io/github/joaquinpelle/Skylight.jl?branch=main)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## General-relativistic ray tracing and radiative transfer in arbitrary spacetimes

Documentation is under construction, and is available [here](https://joaquinpelle.github.io/Skylight.jl/dev).

Skylight works with any spacetime geometry, without the constraints of specific symmetries or the assumption of asymptotic flatness. It is designed with the following goals in mind:
- Fast computational speed
- High accuracy
- Easy extensibility to user-defined spacetimes and radiative models

It uses [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) with [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) to compute the Christoffel symbols from the spacetime metric, and has built-in multithreading parallelism. The integration of the equations is performed with [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). 

For a quick start guide, see [Getting started](https://joaquinpelle.github.io/Skylight.jl/dev/gettingstarted/) (under construction).

### Features

#### Supported spacetimes

* Minkowski spacetime
* Schwarzschild spacetime
* Kerr spacetime
* Johannsen spacetime
* f(R)-Kerr spacetime
* Ruffini-Argüelles-Rueda spacetime for fermionic dark matter
* Boson star spacetimes with quartic self-interaction and solitonic potentials  
* Extensibility to user-defined spacetimes

#### Radiative models

* Shakura-Sunyaev accretion disks
* Geometrically-thin optically-thick accretion disks with user-provided tabulated temperatures
* Ion tori with synchrotron and bremsstrahlung emission 
* Line emission from accretion disks with user-provided emissivity profiles
* Lamppost corona emission and accretion disk illumination profiles
* Extensibility to user-defined radiative models

#### Geometric and dynamical tools

* Spacetime metrics, inverse metrics, volume elements, Christoffel symbols, etc.
* Four-vector scalar products, index raising/lowering, orthogonal projection, normalization, etc.
* Constants of motion in spacetimes with symmetries
* Characteristic radii, like event horizons, ISCOs, etc. 
* Spacetime geodesics integration

#### Radiative transfer mechanisms

* Radiative transfer in vacuum and in emissive/absorptive media

#### Observable quantities

* Bolometric and specific intensities
* Fluxes through arbitrarily oriented surface elements
* Images and spectra
* Generic observation frames (any position and four-velocity)

#### Utilities

* Loading/saving data and configurations from/to HDF5 files
* Units and dimensions management

## Installation

### Requirements
* Julia version at least 1.6

The package is not yet available in the Julia registries. To install it, follow these steps:

1. Clone the repository: `git clone https://github.com/joaquinpelle/Skylight.jl.git`
2. Open the Julia REPL and enter package mode by typing `]`.
3. Add Skylight to your Pkg environment: `dev \path\to\the\repository`
4. Exit package mode with `Ctrl + C` and import Skylight: `using Skylight`

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
