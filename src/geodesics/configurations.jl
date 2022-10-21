export GeodesicsConfigurations

abstract type OutputType end

struct SaveGeodesic <: OutputType end
struct SaveEndpoint <: OutputType end

output_func(output_type::SaveGeodesics) = (sol, i) -> (sol, false)
output_func(output_type::SaveEndpoint) = (sol, i) -> (sol[end], false)

@with_kw struct GeodesicsConfigurations{T<:InitialDataConfigurations, O<:IntegrationOutput, M}
    
    Rmax::Float64
    output::O = SaveEndpoint()
    abstol::Float64 = 1e-21
    reltol::Float64 = 1e-13
    method::M = VCABM()

end

