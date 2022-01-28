function ei(::Type{T}, i::Integer, n::Integer)::Vector{T} where T
    v = zeros(T, n)
    v[i] = one(T)
    return v
end
ei(i::Integer, n::Integer)::Vector{Float64} = ei(Float64, i, n) 


function pivot!(M::Array{T,2}, i::Int, j::Int, Binv::Array{T, 2} = Array{T, 2}(undef, 0, 0); freelines::BitVector = trues(size(M, 1) - 1))::Matrix{T} where T
    updBinv = !iszero(length(Binv))
    m,n = size(M)
    @assert M[i, j] != 0 "pivot  illegal"
    vMi = @view M[i, :]
    updBinv && (vBinv = @view Binv[i, :])
    
    M[i, :] = vMi / M[i, j]
    updBinv && (Binv[i, :] = vBinv/M[i, j])
    
    linenottodo = findall(iszero, @view M[:, j])
    guilylines = findall(!, freelines)
    linenottodo[:] = linenottodo
    for (iter, k) in enumerate(setdiff(1:m, i, linenottodo, guilylines))
        #@show M
        #@show Binv
        M[k, :] -= M[k, j] * vMi
        (k != m) && updBinv && (Binv[iter, :] -= M[k, j] * vBinv)
    end
    return M
end

abignumber(::Type{Float64}) = reinterpret(Float64, (reinterpret(UInt64, typemax(Float64)) - UInt64(1)))/100
abignumber(::Type{Float32}) = reinterpret(Float32, (reinterpret(UInt32, typemax(Float32)) - UInt32(1)))/100
#abignumber(::Type{Float16}) = reinterpret(Float16, (reinterpret(UInt16, typemax(Float16)) - UInt16(1)))/100
abignumber(::Type{Rational{Int128}}) = div(typemax(Int128), Int128(1000))//Int128(1)
abignumber(::Type{Rational{Int64}}) = div(typemax(Int64), Int64(1000))//Int64(1)
abignumber(::Type{Rational{Int32}}) = div(typemax(Int32), Int32(1000))//Int32(1)

tol(::Type{T}) where T <: AbstractFloat = 10*eps(T)
tol(::Type{T}) where T <: Rational = T(0)

function isacanonicalcol(v::AbstractVector)::Bool
    m = length(v)
    return (sum(iszero, v) == (m - 1)) && (sum(isone, v) == 1)
end
function iscanonical(A::AbstractMatrix)::Tuple{Bool, Vector{Int}}
    m, n = size(A)
    basis = zeros(Int, m)
    for i in 1:n
        Ai = @view A[:, i]
        if isacanonicalcol(Ai)
            idx = findfirst(Bool.(Ai))
            basis[idx] = i
        end
    end
    return all(!iszero, basis), basis
end