abstract type AbstractCallbackParameters end

struct Run{C,CB}
    output_data::Array{Float64, 2}
    callback::C
    callback_parameters::CB
    kwargs::Dict
end

@with_kw mutable struct VacuumThreadCache{T}

    point::Vector{Float64} = zeros(4)
    velocity::Vector{Float64} = zeros(4)
    acceleration::Vector{Float64} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T

end

mutable struct VacuumCache{S, C, T}
    
    spacetime::S
    cb_params::C
    multi_thread::Vector{VacuumThreadCache{T}}

end

@with_kw mutable struct NonVacuumThreadCache{T}

    point::Vector{Float64} = zeros(4)
    velocity::Vector{Float64} = zeros(4)
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
    cb_params::C
    τmax::Float64
    observed_energies::Vector{Float64}
    NE::Int
    multi_thread::Vector{NonVacuumThreadCache{T}}

end