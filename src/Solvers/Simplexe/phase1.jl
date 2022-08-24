function phase1!(ss::StandardSimplexe{T}; verbose::Bool = false, timelimit::Float64 = Inf) where  T
    verbose && println("In phase1")
    M = ss.M
    m, n = size(M)
    m -= 1
    n -= 1
    iscanon, basis = iscanonical(@view M[1:end-1, :])
    iscanon && @warn "cannonycal problem feed to phase1!" 
    iscanon && (ss.b_idx[:] = basis)

    Aphase1initial = @view M[1:end-1, 1:end-1]
    base = -ones(Int, m)
    for k in 1:n
        vAk = @view Aphase1initial[:, k]
        iscanoncol = isacanonicalcol(vAk)
        if iscanoncol
            idx = findfirst(t -> abs(t - 1) < 1e-6, vAk)
            base[idx] = k 
        end
    end
    artificialvar = findall(t->t==-1, base)
    @assert length(artificialvar) != 0 "canonycal matrix should have been catched earlier in phase 1 algotithm"
    artificialMatrix = hcat([ei(T, k, m) for k in artificialvar]...)
    base[findall(t->t==-1, base)] = n+1:n+length(artificialvar)
    Aphase1 = [Aphase1initial artificialMatrix]
    cphase1 = [zeros(T, n) ; ones(T, length(artificialvar))]
    bphase1 = M[1:end-1, end]
    ssphase1 = StandardSimplexe(Aphase1, bphase1, cphase1)
    ssphase1.b_idx[:] = base
    for (i, j) in enumerate(base)
        !iszero(ssphase1.M[end, j]) && pivot!(ssphase1.M, i, j)
    end
    verbose && println("simplexe array of phase 1")
    verbose && println(ssphase1)
    status = phase2!(ssphase1, verbose = verbose, timelimit = timelimit)
    verbose && println("solved simplexe array of phase 1")
    verbose && println(ssphase1)
    #if the simplexe array  is optimal but the objective value is different then 0
    if (status == MOI.OPTIMAL) && (ssphase1.vstar > tol(T))
        verbose && println("probleme is unfeasible")
        return MOI.INFEASIBLE
    else
        freelines = (t -> t in 1:n).(ssphase1.b_idx)
        ss.Binv = ssphase1.Binv[freelines, freelines]
        ss.freelines[:] = freelines
        ss.M[1:m, 1:n] = @view ssphase1.M[1:m, 1:n]
        ss.M[1:m, end] = @view ssphase1.M[1:m, end]
        ss.b_idx[:] .= -1
        ss.b_idx[freelines] = ssphase1.b_idx[freelines]
        
        if maximum(base) < n
            println("////////////////////////////////////////////////////////////////////////////////////////////////////////////")
            println("error predicted in base")
            println(ss)
            println("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        end


        for (i,j) in enumerate(base)
            freelines[i] && ss.M[end, j] != zero(T) && pivot!(M, i, j)
        end
        ss.vstar = - ss.M[end, end]
    end
    return MOI.OPTIMAL
end