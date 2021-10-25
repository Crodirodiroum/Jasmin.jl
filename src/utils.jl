ei(::Type{T}, i::Integer, n::Integer) where T = (v = zeros(T, n); v[i] = one(T); v)
ei(i::Integer, n::Integer) = ei(Float64, i, n) 

function pivot!(M::Array{T,2}, i::Int, j::Int) where T
    m,n = size(M)
    @assert M[i, j] != 0 "pivot  illegal"
    M[i, :] = M[i, :] / M[i, j]
    for k in setdiff(1:m, i)
        M[k, :] -= M[k, j] * M[i, :]
    end
    return M
end
