abstract type AbstractStandard{T} end

function Simplexe(A::Array{T, 2}, b::Array{S, 1}, c::Array{W, 1}; method::Type{METHODE} = StandardSimplexe,
        verbose::Bool = false) where {METHODE <: AbstractStandard, T, S, W}
    nT = promote_type(T, S, W)
    !(supertype(nT) <: Real) && (nT = Float64)
    return method(Array{nT, 2}(A), Array{nT, 1}(b), Array{nT, 1}(c), verbose = verbose)
end
    
    
    
include("BaseEnumeration.jl")
include("Simplexe.jl")
