include("BaseEnumeration/main.jl")
include("Simplexe/main.jl")
include("IntervalSimplex/main.jl")
global basesolver = BaseEnumeration
function changeSolver!(::Type{SOLVER} = BaseEnumeration) where {SOLVER <: AbstractSolver}
    global basesolver = SOLVER
end 
function solve!(lp::AbstractLP; verbose::Bool = !lp.issilent)
    solver = basesolver(lp)
    solver(lp, verbose = verbose)
    push!(lp.solvers, solver)
    return lp.status
end