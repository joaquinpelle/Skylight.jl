abstract type AbstractCallbackParameters end

struct Run{C,CB}
    output_data::Matrix{Float64}
    callback::C
    callback_parameters::CB
    kwargs::Dict
end

@with_kw mutable struct VacuumThreadCache{T}
    acceleration::Vector{Float64} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T
end

mutable struct VacuumCache{S, C, T}
    spacetime::S
    cbp::C
    multi_thread::Vector{VacuumThreadCache{T}}
end

@with_kw mutable struct NonVacuumThreadCache{T}
    acceleration::Vector{Float64} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T
    ε::Vector{Float64}
    αε::Vector{Float64}
    jε::Vector{Float64}
end

mutable struct NonVacuumCache{S, M, C, T}
    spacetime::S
    model::M
    cbp::C
    τmax::Float64
    observation_energies::Vector{Float64}
    NE::Int
    multi_thread::Vector{NonVacuumThreadCache{T}}
end