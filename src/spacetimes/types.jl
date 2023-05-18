abstract type AbstractSpacetime end

abstract type AbstractCoordinateSystemClass end

struct CartesianClass <: AbstractCoordinateSystemClass end
struct SphericalClass <: AbstractCoordinateSystemClass end

abstract type ChristoffelCache end