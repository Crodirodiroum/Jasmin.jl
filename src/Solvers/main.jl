include("BaseEnumeration.jl")
include("Simplexe/main.jl")
global basesolver = BaseEnumeration
function changeSolver!(::Type{T}) where T <: AbstractSolver
    global basesolver = T
end 
function solve!(lp::AbstractLP; verbose::Bool = !lp.issilent)
    solver = basesolver(lp)
    solver(lp, verbose = verbose)
    push!(lp.solvers, solver)
    return lp.status
end