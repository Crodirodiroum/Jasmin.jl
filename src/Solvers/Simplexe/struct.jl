mutable struct StandardSimplexe{T, BType <: AbstractMatrix{T}} <: AbstractSolver{T}
    M::Array{T, 2}
    Binv::BType
    xstar::Array{T, 1}
    vstar::T
    b_idx::Array{Int, 1}
    status::Status
    function StandardSimplexe{T, BINV}(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T, BINV}
        #println("in SS")
        @assert !(typeof(T) <: Integer)  "Type $T cannot be a subtype of Integer"
        m, n = size(A)
        M = [A b; c' zero(T)]
        iscanon, b_idx = iscanonical(M)
        if (BINV <: SubArray{T}) 
            if iscanon && all(ss.M[1:end-1, end] .>= 0)
                Binv = @view M[1:end-1, copy(b_idx)]
                b_idx = copy(b_idx)
            else
                M = [A I b; c' zeros(T, 1, m) zero(T)]
                Binv = @view M[1:end-1, n+1:n+m]
            end
        else
            Binv = Array{T, 2}(I, m, m)
            b_idx = -ones(T, m)
        end
        ss = new{T, BINV}()
        ss.M = M
        ss.Binv = Binv
        ss.b_idx = b_idx
        ss.xstar = zeros(T, m)
        ss.status = Unknown
        return ss
    end
end

function StandardSimplexe{T}(A::Matrix{T}, b::Vector{T}, c::Vector{T}) where T
    StandardSimplexe{T, Array{T, 2}}(A, b, c)
end

function StandardSimplexe(lp::AbstractLP{T}) where T
    A = getA(lp)
    b = getb(lp)
    c = getc(lp)
    return StandardSimplexe{T}(A, b, c)
end

function (ss::StandardSimplexe{T})(lp::AbstractLP{T}; verbose::Bool = false) where T
    solve!(ss)
end