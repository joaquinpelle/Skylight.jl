abstract type AbstractCallbackParameters end

struct Run{C,CB}
    output_data::Matrix{Float64}
    callback::C
    callback_parameters::CB
    kwargs::Dict
end

@with_kw mutable struct VacuumThreadCache{SC, CC}
    acceleration::Vector{Float64} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    spacetime_cache::SC
    christoffel_cache::CC
end

mutable struct VacuumCache{S, CBP, SC, CC}
    spacetime::S
    cbp::CBP
    multi_thread::Vector{VacuumThreadCache{SC, CC}}
end

@with_kw mutable struct NonVacuumThreadCache{SC, CC}
    acceleration::Vector{Float64} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    spacetime_cache::SC
    christoffel_cache::CC
    ε::Vector{Float64}
    αε::Vector{Float64}
    jε::Vector{Float64}
end

mutable struct NonVacuumCache{S, M, C, SC, CC}
    spacetime::S
    model::M
    cbp::C
    τmax::Float64
    observation_energies::Vector{Float64}
    NE::Int
    multi_thread::Vector{NonVacuumThreadCache{SC, CC}}
end