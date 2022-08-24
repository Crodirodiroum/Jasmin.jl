mutable struct LinearProblem{T} <: AbstractLP{T}
    A::Matrix{T}
    b::Vector{T}
    c::Vector{T}
    status::MOI.TerminationStatusCode
    solvers::Array{AbstractSolver{T}}
    isMin::Bool
    xstar::Vector{T}
    vstar::T
    issilent::Bool
    timelimit::Float64
    function LinearProblem(A::Matrix{T}, b::Vector{T}, c::Vector{T}; isMin::Bool = true) where T
        isaninteger = (T <: Integer)
        isaninteger && @warn "$T cannot be used by solver, $Float64 used instead"
        timelimit = 60.0
        return new{isaninteger ? Float64 : T}(A, b, c, MOI.OPTIMIZE_NOT_CALLED, [], isMin, zeros(T, size(A, 2)), zero(T), true, timelimit)
    end
    function LinearProblem{T}(A::Matrix, b::Vector, c::Vector; isMin::Bool = true) where T
        isaninteger = (T <: Integer)
        isaninteger && @warn "$T cannot be used by solver, $Float64 used instead"
        timelimit = 60.0
        return new{isaninteger ? Float64 : T}(Array{T, 2}(A), Array{T, 1}(b), Array{T, 1}(c), MOI.OPTIMIZE_NOT_CALLED, [], isMin, zeros(T, size(A, 2)), zero(T), true, timelimit)
    end
end
function LinearProblem(A::Matrix{TA}, b::Vector{Tb}, c::Vector{Tc}; isMin::Bool = true) where {TA, Tb, Tc}
    T = promote_type(TA, Tb, Tc)
    (T <: Integer) && (T = Float64)
    return LinearProblem{T}(A, b, c, isMin = isMin)
end

const LP{T} = LinearProblem{T}

function xstar(lp::LP{T})::Vector{T} where T
    return copy(lp.xstar)
end
function vstar(lp::LP{T})::T where T
    return lp.vstar
end

function nvar(lp::LP)::Int
    return size(lp.A, 2)
end
function ncons(lp::LP)::Int
    return size(lp.A, 1)
end
function getA(lp::LP{T})::Array{T, 2} where T
    return lp.A
end
function getb(lp::LP{T})::Array{T, 1} where T
    return lp.b
end
function getc(lp::LP{T})::Array{T, 1} where T
    return lp.c
end
#=
function addconstraint!(lp::LinearProblem{T}, v::Vector{T}, b::T) where T
    lp.A = [lp.A; v']
    lp.b = [lp.b; b]
    #lp.c::Vector{T}
    lp.status = Unknown
    addconstraint!.(lp.solvers, v, b)
    lp.xstar[:] .= zero(T)
    lp.vstar = lp.isMin ? abignumber(T) : - abignumber(T)
end
function addvariable!(lp::LinearProblem{T}) where T
    lp.A = [lp.A zeros(T, size(lp.A, 1))]
    lp.c = [lp.c; zero(T)]
    addvariable!.(lp.solvers)
end
=#