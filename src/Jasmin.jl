module Jasmin

using LinearAlgebra, Combinatorics, Test

import Base.println

#=
export AbstractStandard, BaseEnumeration, StandardSimplexe, solve!,  Simplexe,
     AbstractStatus, Optimal, Infeasible, Unknown, Unbounded,
     xstar, vstar
=#
include("utils.jl")
include("Status.jl")

abstract type AbstractSolver{T} end
abstract type AbstractLP{T} end

include("LinearProblemRepresentation/main.jl")
include("Solvers/main.jl")

#include("test/main.jl")


#by default, Jasmin does base enumeration to solve a problem. 
changeSolver!(basesolver)
end # module
