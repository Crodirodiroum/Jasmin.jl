mutable struct SimplexeDictionary{T} <: AbstractSolver{T}
    xstar::Array{T, 2}
    vstar::T
    colindexs::Array{Int, 2}
    removedLines::BitVector
    n::Int
    function SimplexeDictionary(lp::AbstractLP{T}) where T
        n = nvar(lp)
        m = ncons(lp)
        xstar = Array{T, 2}(undef, n + 1, 0) #a made up variable
        vstar = typemax(T)
        colindexs = Array{Int, 2}(undef, m + 1, 0) #made up constraint
        removedLines = trues(m)
        return new{T}(xstar, vstar, colindexs, removedLines, n)
    end
end
const SD{T} = SimplexeDictionary{T}
function (be::SD{T})(lp::AbstractLP{T}; verbose::Bool = !lp.issilent) where T
    A = getA(lp)
    b = getb(lp)
    c = getc(lp)


    #todo
    lp.status = MOI.Infeasible
    lp.xstar = zeros(T,0)
    lp.vstar = zero(T)
    return lp.status
end