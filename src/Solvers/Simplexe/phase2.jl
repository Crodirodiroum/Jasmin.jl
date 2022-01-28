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
        if (ss.freelines[k]) && (M[k, entering] > 0) && ((tobeat = b[k]/M[k, entering]) < val)
            leaving = k
            val = tobeat
        end
    end
    verbose && @show leaving
    return leaving, entering
end

function phase2!(ss::StandardSimplexe{T}; verbose::Bool = false, kmax = 10000) where T
    k = 0
    m2, n2 = size(ss.M)
    m, n = m2 - 1,  n2 -1
    while !isOptimal(ss) && (k+=1) <= kmax
        verbose && @show k
        verbose && println(ss)
        i, j = findpivot(ss, verbose = verbose)
        if i == -1
            ss.status = Unbounded
            break
        end
        pivot!(ss.M, i, j, ss.Binv, freelines = ss.freelines)
        ss.b_idx[i] = j
    end
    if isOptimal(ss) 
        ss.status = Optimal
        ss.xstar = zeros(T, n)
        freebasis = ss.b_idx[ss.freelines]
        xB = ss.M[[ss.freelines; false], end]
        ss.xstar[freebasis] = xB
        ss.vstar = -ss.M[end, end]
    end  
    
    status = ss.status
    verbose && (println("$status Simplexe");  println(ss))
    verbose && @show ss.xstar
    verbose && @show ss.vstar
    ss.status
end