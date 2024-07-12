# Configurations

In Skylight.jl, the configurations objects reunite the main components of the astrophysical problems: the spacetime, the radiative model, and the camera. Additionally, it contains a parameter to completely determine the geometrized code units together with the $G = c = 1$ conditions. There are different configuration types depending on the method (observer-to-emitter or emitter-to-observer) and whether the radiative transport occurs in a vacuum, which allows for certain simplifications. However, the emitter-to-observer method is not fully implemented yet in this version of Skylight.jl, so we only consider the observer-to-emitter method. 

## Examples

For vacuum transport problems using the observer-to-emitter method, the configurations are implemented by the [`VacuumOTEConfigurations`](@ref) type. They can be constructed as follows:

```julia
configurations = VacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0)
```

where `spacetime` is the spacetime object, `camera` is the camera object, and `model` is the radiative model object. The `unit_mass_in_solar_masses` is equal to $M/M_\odot$ where $M=1$ is the unit mass of the geometrized unit system of the code. An important point is that this parameter must equal the value set in any radiative model object that requires it, as the [`NovikovThorneDisk`](@ref) model for example. We acknowledge this is a non-ideal solution, and we are working on a more consistent way to handle the units in Skylight.jl. The advantage of this method is that it avoids re-running ray-tracing when the absolute unit scale changes.

For non-vacuum transport problems, such as those involving an [`IonTorus`](@ref) radiative model, the configurations are implemented by the [`NonVacuumOTEConfigurations`](@ref) type. The can be constructed as follows:

```julia
configurations = NonVacuumOTEConfigurations(spacetime = spacetime,
    camera = camera,
    radiative_model = model,
    unit_mass_in_solar_masses = 1.0,
    observation_energies = [1e-8, 1e-7])
```

The `observation_energies` parameter specifies the (CGS) energies at which the specific intensity is to be computed in the observation frame. In non-vacuum problems, this has to be specified right away, as the solution depends on the energy in a non-trivial manner, as opposed to vacuum problems where only the invariant $I_\nu / \nu^3$ needs to be connected between the observation and emission frames.

## Types

```@docs
VacuumOTEConfigurations
NonVacuumOTEConfigurations
```