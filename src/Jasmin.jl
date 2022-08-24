module Jasmin

using LinearAlgebra, Combinatorics, Test, MathOptInterface
const MOI = MathOptInterface


#=
export AbstractStandard, BaseEnumeration, StandardSimplexe, solve!,  Simplexe,
     xstar, vstar
=#
include("utils.jl")

abstract type AbstractSolver{T} end
abstract type AbstractLP{T} end

include("LinearProblemRepresentation/main.jl")
include("Solvers/main.jl")
include("Wrapper/main.jl")
#include("test/main.jl")


#by default, Jasmin does base enumeration to solve a problem. 
changeSolver!()
end # module
