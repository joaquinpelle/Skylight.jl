# Skylight
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joaquinpelle.github.io/Skylight.jl/dev)
[![Build Status](https://github.com/joaquinpelle/Skylight.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/joaquinpelle/Skylight.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/joaquinpelle/Skylight.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/joaquinpelle/Skylight.jl)
[![Coverage](https://coveralls.io/repos/github/joaquinpelle/Skylight.jl/badge.svg?branch=main)](https://coveralls.io/github/joaquinpelle/Skylight.jl?branch=main)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A general-relativistic ray tracing and radiative transport code for arbitrary spacetimes. With multithreading and automatic differentiation support.

> [Reference](https://academic.oup.com/mnras/article-abstract/515/1/1316/6631564): Pelle, Joaquin, Oscar Reula, Federico Carrasco, and Carlos Bederian. "Skylight: a new code for general-relativistic ray-tracing and radiative transfer in arbitrary space–times." Monthly Notices of the Royal Astronomical Society 515, no. 1 (2022): 1316-1327.

### Requirements

At least Julia v1.6

### Folder contents

    .
    ├── docs                    # Documentation files
    │   ├── TOC.md              # Table of contents
    │   ├── faq.md              # Frequently asked questions
    │   ├── misc.md             # Miscellaneous information
    │   ├── usage.md            # Getting started guide 
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

If you use this software in your work, please cite it using the following metadata

```
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