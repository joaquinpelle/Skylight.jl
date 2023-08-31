# Catalogue of spacetimes

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractSpacetime && isconcretetype(t)
```

```@docs
Skylight.RARSpacetime
```