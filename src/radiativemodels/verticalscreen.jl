#Does not need to be documented
@with_kw struct VerticalScreen <: AbstractSurfaceEmissionModel
    x::Float64
    horizontal_side::Float64
    vertical_side::Float64
end

stationarity(::VerticalScreen) = IsStationary()