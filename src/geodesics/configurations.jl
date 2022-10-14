export GeodesicsConfigurations

@with_kw struct GeodesicsConfigurations{T<:InitialDataConfigurations,M}
    
    initial_data_configs::T
    abstol::Float64
    reltol::Float64
    method::M
    Rmax::Float64

end
