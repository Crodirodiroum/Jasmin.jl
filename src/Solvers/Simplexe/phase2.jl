function findpivot(ss::StandardSimplexe{T}; verbose::Bool = false)::Tuple{Int, Int} where T
    if isdegenerated(ss)
        return BlandFindpivot(ss, verbose = verbose)
    end
    M = ss.M
    m = size(M, 1) - 1
    c = @view M[end, 1:end-1]
    b = @view M[1:end-1, end]
    entering = argmin(c)
    verbose && @show entering
    @assert c[entering] < 0 "simplexe is in an optimal state" #ss was optimal
    val = T(Inf)
    leaving = -1
    for k in 1:m
        if (ss.freelines[k]) && (M[k, entering] > 0) && ((tobeat = b[k]/M[k, entering]) < val)
            leaving = k
            val = tobeat
        end
    end
    verbose && @show leaving
    return leaving, entering
end
function BlandFindpivot(ss::StandardSimplexe{T}; verbose::Bool = verbose)::Tuple{Int, Int} where T
    M = ss.M
    m = size(M, 1) - 1
    c = @view M[end, 1:end-1]
    b = @view M[1:end-1, end]
    entering = findfirst(t->t<zero(T), c)
    verbose && @show entering
    @assert c[entering] < 0 "simplexe is in an optimal state" #ss was optimal
    val = T(Inf)
    leaving = -1
    for k in 1:m
        if (ss.freelines[k]) && (M[k, entering] > 0) && ((tobeat = b[k]/M[k, entering]) < val)
            leaving = k
            val = tobeat
        end
    end
    verbose && @show leaving
    return leaving, entering
end
function isdegenerated(ss::StandardSimplexe)::Bool
    M = ss.M
    vTD = @view M[1:end-1, end]
    return any(iszero, vTD[ss.freelines])
end

function phase2!(ss::StandardSimplexe{T}; verbose::Bool = false, kmax = 10000, timelimit::Float64 = Inf)::MOI.TerminationStatusCode where T
    timephase2start = time()
    status = MOI.OTHER_ERROR
    k = 0
    m2, n2 = size(ss.M)
    m, n = m2 - 1,  n2 -1
    while !isOptimal(ss) && (k+=1) <= kmax
        (time() - timephase2start > timelimit) && return MOI.TIME_LIMIT
        verbose && @show k
        verbose && println(ss)
        i, j = findpivot(ss, verbose = verbose)
        if i == -1
            return MOI.DUAL_INFEASIBLE
        end
        pivot!(ss.M, i, j, ss.Binv, freelines = ss.freelines)
        ss.b_idx[i] = j
    end
    if isOptimal(ss) 
        status = MOI.OPTIMAL
        ss.xstar = zeros(T, n)
        freebasis = ss.b_idx[ss.freelines]
        xB = ss.M[[ss.freelines; false], end]
        ss.xstar[freebasis] = xB
        ss.vstar = -ss.M[end, end]
    end  
    
    verbose && (println("$status Simplexe");  println(ss))
    verbose && @show ss.xstar
    verbose && @show ss.vstar
    status
end