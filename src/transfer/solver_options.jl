export SolverOptions
export SaveEndstate, SaveSolution

abstract type OutputType end

struct SaveSolution <: OutputType end
struct SaveEndstate <: OutputType end

output_func(output_type::SaveSolution) = (sol, i) -> (sol, false)
output_func(output_type::SaveEndstate) = (sol, i) -> (sol[end], false)

@with_kw struct SolverOptions{M, O}
    
    method::M = VCABM()
    output_type::O = SaveEndstate()
    abstol::Float64 = 1e-21
    reltol::Float64 = 1e-13

end