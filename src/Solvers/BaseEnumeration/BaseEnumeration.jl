mutable struct BaseEnumeration{T} <: AbstractSolver{T}
    xstar::Array{T, 2}
    vstar::T
    colindexs::Array{Int, 2}
    removedLines::BitVector
    n::Int
    function BaseEnumeration(lp::AbstractLP{T}) where T
        n = nvar(lp)
        m = ncons(lp)
        xstar = Array{T, 2}(undef, n + 1, 0) #a made up variable
        vstar = typemax(T)
        colindexs = Array{Int, 2}(undef, m + 1, 0) #made up constraint
        removedLines = trues(m)
        return new{T}(xstar, vstar, colindexs, removedLines, n)
    end
end
const BE{T} = BaseEnumeration{T}
function (be::BE{T})(lp::AbstractLP{T}; verbose::Bool = !lp.issilent) where T
    time0 = time()
    A = getA(lp)
    b = getb(lp)
    c = getc(lp)
    r1 = rank(A)
    r2 = rank([A b])
    if r1 != r2 
        verbose && println("Infeasible since rank(A) != rank([A b])")
        #this is Infeasible
        lp.status = MOI.Infeasible
    end
    if r1 < size(A, 1)
        #there are redundent constraint
        removedLines = reduceproblem(A, b)
        @assert sum(removedLines) == r1 "Impossible to reduce the problem to a feasible form"
        be.removedLines[:] = removedLines
        verbose && println("A not full rank, $(sum(!, be.removedLines)) rows removed from problem.")
        A = A[removedLines, :]
        b = b[removedLines]
    end
    m, n = size(A)
    A = [A zeros(T, m); c' -one(T)]
    b = [b ; -abignumber(T)]
    c = [c; zero(T)]
    verbose && @show A
    verbose && @show b
    verbose && @show c
    m, n = size(A)
    C =  combinations(1:n, m)
    verbose && println("will loop over the $(length(C)) possible basis.")
    for (iter, comb) in enumerate(C)
        if time() - time0 > lp.timelimit
            lp.status = MOI.TIME_LIMIT
            return MOI.TIME_LIMIT
        end
        B = @view A[:, comb]
        if rank(B) == m
            xB = B\b
            if all(xB .>= - tol(T))
                cB = c[comb]
                vtemp = dot(cB, xB)
                if vtemp < be.vstar - tol(T)
                    #found a better solution
                    be.vstar = vtemp
                    be.xstar = zeros(T, be.n + 1, 1)
                    be.xstar[comb, end] = xB
                    be.colindexs = reshape(comb, m, 1)
                elseif be.vstar - tol(T) <= vtemp <= be.vstar + tol(T)
                    #found a similar solution
                    be.xstar = [be.xstar zeros(T, be.n + 1)]
                    be.xstar[comb, end] = xB
                    be.colindexs = [be.colindexs comb]
                end
                verbose && println("basis $iter feasible, solution $vtemp")
            else
                verbose && println("basis $iter full rank but not feasible xB = $xB")
            end
        else
            verbose && println("basis $iter not full rank")
        end
    end
    if be.vstar == typemax(T)
        lp.status = MOI.INFEASIBLE
    elseif be.vstar == -abignumber(T)
        lp.status = MOI.DUAL_INFEASIBLE
    else
        lp.status = MOI.OPTIMAL
        lp.xstar = be.xstar[1:end - 1, 1]
        lp.vstar = be.vstar
    end
    return lp.status
end

function reduceproblem(A::Matrix{T}, b::Vector{T})::BitVector where T
    M = [A b]
    m, n = size(M)
    reducedprob = trues(m)
    for i in 1:m
        j = findfirst(!iszero, M[i, :])
        @assert j != n "unexpected error in reduceproblem, line $i, M = $M"
        isnothing(j) ? (reducedprob[i] = false) : pivot!(M, i, j)
    end
    return reducedprob
end
function getdual(lp::AbstractLP{T}, be::BaseEnumeration{T}) where T
    baseindex = be.colindexs[1:end - 1, 1]
    removedlines = be.removedLines[1:end]
    A = getA(lp)
    b = getb(lp)
    c = getc(lp)
    m, n = size(A)

    B = @view A[removedlines, baseindex]
    cB = @view c[baseindex]

    lambda = zeros(T, m)
    lambda[removedlines] = (B')\cB
    return lambda
end