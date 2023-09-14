# Catalogue of spacetimes

### Minkowski spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractMinkowskiSpacetime && isconcretetype(t)
```

### Schwarzschild spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractSchwarzschildSpacetime && isconcretetype(t)
```

### Kerr spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractKerrSpacetime && isconcretetype(t)
```

### Johannsen spacetime

```@docs
Skylight.JohannsenSpacetime
```

### f(R)-Kerr spacetime

```@docs
Skylight.FRKerrSpacetime
```

### Charged wormhole spacetime

```@autodocs
Modules = [Skylight]
Filter = t -> typeof(t) === DataType && t <: Skylight.AbstractChargedWormholeSpacetime && isconcretetype(t)
```

### RAR spacetime

```@docs
Skylight.RARSpacetime
```

### Boson star spacetimes

```@docs
Skylight.BosonStarSpacetime
```