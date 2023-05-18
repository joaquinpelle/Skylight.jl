abstract type CallbackParameters end

struct Run{C,CB}
    output_data::Array{Float64, 2}
    callback::C
    callback_parameters::CB
    kwargs::Dict
end

@with_kw mutable struct VacuumThreadCache{T}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T

end

mutable struct VacuumCache{S, C, T}
    
    spacetime::S
    cb_params::C
    multi_thread::Array{VacuumThreadCache{T},1}

end

@with_kw mutable struct NonVacuumThreadCache{T}

    point::Array{Float64, 1} = zeros(4)
    velocity::Array{Float64, 1} = zeros(4)
    acceleration::Array{Float64, 1} = zeros(4)
    christoffel::Array{Float64, 3} = zeros(4,4,4)
    christoffel_cache::T
    ε::Array{Float64, 1}
    αε::Array{Float64, 1}
    jε::Array{Float64, 1}

end

mutable struct NonVacuumCache{S, M, C, T}
    
    spacetime::S
    model::M
    cb_params::C
    τmax::Float64
    observed_energies::Array{Float64, 1}
    NE::Int64
    multi_thread::Array{NonVacuumThreadCache{T},1}

end