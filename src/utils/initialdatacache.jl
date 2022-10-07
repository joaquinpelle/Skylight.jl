abstract type InitialDataCache end

@with_kw mutable struct OTEInitialDataCache
    gμν::Matrix{Float64} = zeros(4,4)
    tμ::Vector{Float64} = zeros(4)
end

@with_kw mutable struct ETOInitialDataCache
    gμν::Matrix{Float64} = zeros(4,4)
    tμ::Vector{Float64} = zeros(4)
    eμa::Matrix{Float64} = zeros(4,4)
end

function dump_vector_in!(cache, vector)
    cache.tμ = vector
end

function dump_∂t_in!(cache)
    cache.tμ = ∂t()
end

function dump_metric_in!(cache, position, spacetime::Spacetime)
    set_metric!(cache.gμν, position, spacetime)
end

function unpack_views(cache::OTEInitialDataCache)

    @views begin
        gμν = cache.gμν
        tμ = cache.tμ
    end

    return gμν, tμ

end

function unpack_views(cache::ETOInitialDataCache)

    @views begin
        gμν = cache.gμν
        tμ = cache.tμ
        eμa = cache.eμa
    end

    return gμν, tμ, eμa

end