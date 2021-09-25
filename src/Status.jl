abstract type AbstractStatus end
struct Optimal <: AbstractStatus end
struct Infeasible <: AbstractStatus end
struct Unknown <: AbstractStatus end
struct Unbounded <: AbstractStatus end