#=
#code to add variables to the model and the given solver
function canaddvariables(SO::Type{StandardSimplexe{T, Matrix{T}}}) where T
    return true
end
function addvariables!(ss::StandardSimplexe{T, Matrix{T}}, k::Int) where T
    ss.M = [ss.M[:, 1:end-1] zeros(T, m, k), ss.M[:, end]]
    ss.xstar = [ss.xstar; zero{T}]
end

#code to add constraints to the model and the given solver
function canaddconstraint(SO::Type{StandardSimplexe{T, Matrix{T}}}) where T
    return false #addconstraint! not finished
end

function addconstraint!(ss::StandardSimplexe{T, Matrix{T}}, ai::Vector{T}, b::T)
    ss.M = [ss.M[1:end-1, :]; ai' bi; ss.M[end, :]]
    error("addconstraint! is not finished")
    if ss.status == Infeasible() #adding constraint wont make the problem feasible
        return ss.status
    end
    iscanon, b_idx = iscanonical(M)
    m = size(ss.Binv, 1)
    if iscanon
        ss.b_idx = [ss.b_idx; b_idx[end]]
        Binv = [ss.Binv zeros(T, m); zeros(T, 1, m) 1]
        ss.xstar[b_idx] = M[1:end-1, end]
        ss.status = genstatus(ss)
        return status
    end
    return Unknown()
end

#code to change constraints to the model and the given solver
function canchangeconstraint(SO::Type{StandardSimplexe{T, Matrix{T}}}) where T
    return true
end
function changeConstraint!(ss::StandardSimplexe, b::Vector{T})
    ss.M[1:end-1, end] = ss.Binv*b
    return genstatus(ss)
end
=#
function findfirstbasis(ss::StandardSimplexe{T}; verbose::Bool = false) where T
    verbose && println("in findfirstbasis")
    m, n = size(ss.A)
    for k in 1:m
        if ss.M[k, end] <  0
            ss.M[k, :] *= -1
        end
    end
    M, b_idx = phase1(copy(ss.A), copy(ss.b), copy(ss.c), verbose = verbose)
    verbose && b_idx == [-1] && println("phase 1 did not find a feasible basis")
    ss.b_idx = b_idx
    ss.M[:, :] = M
    ss.status = (b_idx == [-1]) ? Infeasible : Unknown
end

function phase1(A::Matrix{T},b::Vector,c::Vector; verbose::Bool = false) where  T
    verbose && println("In phase1")
    m, n = size(A)
    Aphase1initial = copy(A)
    bphase1 = b 
    base = -ones(Int, m)
    for k in 1:m 
        if bphase1[k] < 0
            bphase1[k] *= -1
            Aphase1initial[k, :] *= -1
        end
    end
    for k in 1:n
        vAk = @view Aphase1initial[:, k]
        iscanoncol = isacanonicalcol(vAk)
        if isiso
            base[idx] = findfirst(vAk)
        end
    end
    artificialvar = findall(t->t==-1, base)
    artificialMatrix = hcat([ei(T, k, m) for k in artificialvar]...)
    if length(artificialvar) == 0
        M = [Aphase1initial bphase1; c' zero(T)]
        for (k, bk) in enumerate(base)
            pivot!(M, k, bk)
        end
        return M, base
    end
    base[findall(t->t==-1, base)] = n+1:n+length(artificialvar)
    @show base
    Aphase1 = [Aphase1initial artificialMatrix]
    
    
    cphase1 = [zeros(T, n) ; ones(T, length(artificialvar))]
    ssphase1 = StandardSimplexe(Aphase1, bphase1, cphase1, verbose = verbose, b_idx = base)
    for (k, bk) in enumerate(base)
        pivot!(ssphase1.M, k, bk)
    end
    verbose && println("simplexe array of phase 1")
    verbose && println(ssphase1)
    solve!(ssphase1, verbose = verbose)
    verbose && println("solved simplexe array of phase 1")
    verbose && println(ssphase1)
    #if the simplexe array  is optimal but the objective value is different then 0
    if (ssphase1.status == Optimal()) && (ssphase1.vstar != zero(T))
        println("probleme is unfeasible")
        return [ssphase1.M[1:m, 1:n] ssphase1.M[1:m, end]; c' zero(T)], [-1]
    else
        M = [ssphase1.M[1:m, 1:n] ssphase1.M[1:m, end]; c' zero(T)]
        base = ssphase1.b_idx
        @assert issubset(base, 1:n) "artificial variable in base of phase 1"
        for (i,j) in enumerate(base)
            pivot!(M, i, j)
        end

        return M, base
    end
end
                                
function isOptimal(ss::StandardSimplexe{T}) where T
    return all(ss.M[end, 1:end-1] .>= -10*eps())
end                     


function solve!(ss::StandardSimplexe{T}; verbose::Bool = false, kmax = 1000) where T
    ss.b_idx == [-1] && findfirstbasis(ss, verbose = verbose)
    k = 0
    if ss.status == Infeasible()
        return ss.status
    end
    m2, n2 = size(ss.M)
    m, n = m2 - 1,  n2 -1
    while !isOptimal(ss) && (k+=1) <= kmax
        verbose && println(ss)
        i, j = findpivot(ss, verbose = verbose)
        if i == -1
            ss.status = Unbounded()
            break
        end
        pivot!(ss.M, i, j)
        ss.b_idx[i] = j
    end
    if isOptimal(ss) 
        ss.status = Optimal()
        ss.xstar = zeros(T, n)
        ss.xstar[ss.b_idx] = ss.M[1:end-1, end]
        ss.vstar = -ss.M[end, end]
    end  
    status = ss.status
    verbose && (println("$status Simplexe");  println(ss))
    ss.Binv = ss.A[:, ss.b_idx]^(-1)
    ss.status
end

function findpivot(ss::StandardSimplexe{T}; verbose::Bool = false) where T
    M = ss.M
    m = size(M, 1) - 1
    c = ss.M[end, 1:end-1]
    b=  ss.M[1:end-1, end]
    entering = argmin(c)
    verbose && @show entering
    @assert c[entering] < 0 "simplexe is in an optimal state" #ss was optimal
    val = T(Inf)
    leaving = -1
    for k in 1:m
        if (M[k, entering] > 0) && ((tobeat = b[k]/M[k, entering]) < val)
            leaving = k
            val = tobeat
        end
    end
    verbose && @show leaving
    return leaving, entering
end
    
    
function xstar(ss::StandardSimplexe{T}) where T
    ss.status != Optimal && @warn "probleme not optimal"
    return copy(ss.xstar)
end
function vstar(ss::StandardSimplexe{T}) where T
    ss.status != Optimal && @warn "probleme not optimal"
    return ss.vstar
end
    

function println(ss::StandardSimplexe)
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
        println(prd)
    end
end
