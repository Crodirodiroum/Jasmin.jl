mutable struct InteriorPoint{T} <: AbstractSolver{T}
    xstar::Array{T, 2}
    vstar::T
    n::Int
    function InteriorPoint(lp::AbstractLP{T}) where T
        n = nvar(lp)
        m = ncons(lp)
        xstar = Array{T, 2}(undef, n, 0) #a made up variable
        vstar = typemax(T)
        return new{T}(xstar, vstar, n)
    end
end
const IP{T} = InteriorPoint{T}
function (ip::IP{T})(lp::AbstractLP{T}; verbose::Bool = !lp.issilent) where T
    time0 = time()
    A = getA(lp)
    b = getb(lp)
    c = getc(lp)

    #TODO
    #solve min c'x
    #subject to A*x = b
    #x[i] >= 0 foralll i

    be.vstar = zero(T) #TODO
    be.xstar = zeros(T, 0) #TODO
    lp.status = MOI.INFEASIBLE #TODO
    return lp.status
end