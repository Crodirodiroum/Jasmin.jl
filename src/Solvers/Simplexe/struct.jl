mutable struct StandardSimplexe{T} <: AbstractSolver{T}
    M::Array{T, 2}
    Binv::Array{T, 2}
    xstar::Array{T, 1}
    vstar::T
    b_idx::Array{Int, 1}
    freelines::BitVector#index k is true if line k is linearly independent of the previous lines, dependent lines are found during phase 1.
    function StandardSimplexe{T}(M::Matrix{T}) where {T}
        @assert !(typeof(T) <: Integer)  "Type $T cannot be a subtype of Integer"
        m, n = size(M).-1
        ss = new{T}()
        ss.M = M
        ss.Binv = Array{T, 2}(I, m, m)
        ss.b_idx = -ones(T, m)
        ss.xstar = zeros(T, n)
        ss.freelines = trues(m) #lines are inocent until proven guilty
        return ss
    end
end
function StandardSimplexe{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    M = [A b; c' zero(T)]
    return  StandardSimplexe{T}(M)
end
function StandardSimplexe(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    M = [A b; c' zero(T)]
    return  StandardSimplexe{T}(M)
end
function StandardSimplexe(lp::AbstractLP{T}) where T
    A = getA(lp)
    b = getb(lp)
    c = getc(lp)
    return StandardSimplexe{T}(A, b, c)
end

function (ss::StandardSimplexe{T})(lp::AbstractLP{T}; verbose::Bool = !lp.issilent) where T
    verbose = true #todo change back
    verbose && println("got a problem to solve, heres is it's linear Representation \n", lp)
    timeoptimizationstart = time()
    m = size(ss.M, 1) - 1
    for k in 1:m
        if ss.M[k, end] <  0
            ss.Binv[k, :] *= -1
            ss.M[k, :] *= -1
        end
    end
    M = ss.M
    vM = @view M[1:end-1, :]
    iscanon, basis = iscanonical(vM)
    status = MOI.OPTIMIZE_NOT_CALLED
    if iscanon
        ss.b_idx[:] = basis
        isOptimal(ss) && (lp.status = MOI.OPTIMAL; lp.xstar[:] = ss.xstar; ss.vstar = lp.vstar; return lp.status)
        status = phase2!(ss, verbose = verbose, timelimit = lp.timelimit)
    else
        statusphase1 = phase1!(ss, verbose = verbose, timelimit = lp.timelimit)
        (statusphase1 == MOI.INFEASIBLE) && (lp.status = MOI.INFEASIBLE; return lp.status) 
        isOptimal(ss) && (lp.status = MOI.OPTIMAL; lp.xstar[:] = ss.xstar; ss.vstar = lp.vstar; return lp.status)
        currenttime = time()
        timelimit = lp.timelimit - (currenttime - timeoptimizationstart)
        (statusphase1 == MOI.OPTIMAL) && (status = phase2!(ss, verbose = verbose, timelimit = timelimit))
    end
    lp.status = status
    isOptimal(ss) && (lp.xstar[:] = ss.xstar; lp.vstar = ss.vstar)
    verbose && @show lp.vstar
    verbose && @show lp.xstar
    return lp.status
end


function Base.print(ss::StandardSimplexe)
    m, n = size(ss.M)
    prd = "+" * prod("-------+" for _ in 1:n+1)
    println(prd)
    print("| V.B.  |")
    for k in 1:n-1
        s = ""
        xk = " x_$k "
        s*= xk
        s*= prod(" " for _ in 1:7-length(xk))
        print(s)
        print("|")
    end
    #@show ss.b_idx
    print(" T.D.  |")
    println()
    println(prd)
    for i in 1:m
        if i != m
            s = "| "#------"
            xk = ss.b_idx != [-1] ? "x_$(ss.b_idx[i])" : "   "
            s*= xk
            s*= prod(" " for _ in 1:6-length(xk))
            print(s)
        else
            print("| -z    ")
        end
        for j in 1:n
            print("| ")
            str = string(Float64(ss.M[i, j]))
            if length(str) >= 5
                str = str[1:5]
            end
            print(str)
            if length(1:(5 - length(str))) != 0
                print(prod(" "  for _ in 1:(5 - length(str))))
            end
            print(" ")
        end
        println("|")
        print(prd)
    end
end
function xstar(ss::StandardSimplexe{T})::Vector{T} where T
    #ss.status != Optimal && @warn "probleme not optimal"
    return copy(ss.xstar)
end
function vstar(ss::StandardSimplexe{T})::T where T
    #ss.status != Optimal && @warn "probleme not optimal"
    return ss.vstar
end
function isOptimal(ss::StandardSimplexe{T})::Bool where T
    return all(ss.M[end, 1:end-1] .>= -tol(T))
end   