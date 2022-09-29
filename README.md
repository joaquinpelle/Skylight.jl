# Skylight
A general-relativistic ray tracing and radiative transport code for arbitrary spacetimes

> [Reference](https://academic.oup.com/mnras/article-abstract/515/1/1316/6631564): Pelle, Joaquin, Oscar Reula, Federico Carrasco, and Carlos Bederian. "Skylight: a new code for general-relativistic ray-tracing and radiative transfer in arbitrary space–times." Monthly Notices of the Royal Astronomical Society 515, no. 1 (2022): 1316-1327.

### Purpose

This is a minimal version of `skylight`, currently under development for adding the following features: 
  * Interfacing with 3D GRMHD simulations data
  * Spacetimes from full numerical relativity
  
### Requirements

### Folder contents

> Adapted from [this guide](https://github.com/kriasoft/Folder-Structure-Conventions)

    .
    ├── data                    # Data files
    │   ├── io                  # Input and output 
    │   ├── meta                # Metadata files
    ├── docs                    # Documentation files
    │   ├── TOC.md              # Table of contents
    │   ├── faq.md              # Frequently asked questions
    │   ├── misc.md             # Miscellaneous information
    │   ├── usage.md            # Getting started guide 
    ├── run                     # Scripts and notebooks to run the code
    ├── src                     # Source files
    ├── test                    # Test files 
    │   ├── benchmarks          # Load and stress tests
    │   ├── integration         # End-to-end, integration tests
    │   └── unit                # Unit tests


