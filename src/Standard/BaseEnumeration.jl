mutable struct BaseEnumeration{T} <: AbstractStandard{T}
    A::Array{T, 2}
    b::Array{T, 1}
    c::Array{T, 1}
    xstar::Array{T, 1}
    vstar::T
    status::AbstractStatus
    function BaseEnumeration(A::Array{T, 2}, b::Array{T, 1}, c::Array{T, 1}) where T
        m,n = size(A)
        @assert !(typeof(T) <: Integer)  "Type $T cannot be a subtype of Integer"
        @assert length(b) == m "dimension of A and b mismatch, size(A) = ($m, $n), length(b) = $(length(b)) != $m"
        @assert length(c) == n "dimension of A and c mismatch, size(A) = ($m, $n), length(c) = $(length(c)) != $n"
        @assert rank(A) == m "A is not a full rank Matrix"
        be = new{T}(A, b, c, [],  T(Inf), Unknown())
        return be
    end
end

function solve!(be::BaseEnumeration{T}; verbose::Bool = false) where T
    m, n = size(be.A)
    vopt = T(Inf)
    sol = zeros(T, n)
    allComb = [c for c in combinations(1:n, m)]
    for comb in allComb
        verbose && println(comb)
        B = be.A[:, comb]
        if rank(B) == m && all((x_B = B\be.b) .>= 0) && (v_B = dot(be.c[comb], x_B)) < vopt
            vopt = v_B
            sol[:] .= zero(T)
            sol[comb] = x_B
        end
    end
    if vopt ==  T(Inf)
        be.status = Infeasible()
    else
        be.status = Optimal()
        be.vstar = vopt
        be.xstar = sol
    end
    be.status
end
            