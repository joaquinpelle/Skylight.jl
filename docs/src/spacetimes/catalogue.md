# Catalogue of spacetimes

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractSpacetime && isconcretetype(t)
```
### Abstract types
```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractSpacetime && isabstracttype(t)
```

