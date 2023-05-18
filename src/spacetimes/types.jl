abstract type AbstractSpacetime end

abstract type CoordinateSystemClass end

struct CartesianClass <: CoordinateSystemClass end
struct SphericalClass <: CoordinateSystemClass end

abstract type ChristoffelCache end