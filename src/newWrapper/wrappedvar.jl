mutable struct WrappedListofvar{T}
    lower::Array{T, 1}
    upper::Array{T, 1}
    indexs::Array{Tuple{Int, Int}, 1}
    haslower::BitVector
    hasupper::BitVector
    isequalto::BitVector
    function Listofvar{T}()
        return new{T}([], [], [], falses(0), falses(0), falses(0))
    end
end
function numberofvarandconstraint(wlv::WrappedListofvar)::Tuple{Int, Int}
    nvar = 0
    ncons = 0
    for i in 1:length(wlv.lower)
        if wlv.haslower[i] && wlv.lower[i] <= zero(T)
            nvar += 1
        else 
            nvar += 2
        end

        if wlv.haslower[i] && wlv.lower[i] != zero(T)
            ncons += 1
            nvar += 1
        end
        if lv.hasupper[i]
            ncons += 1
            nvar += 1
        end
    end
    return nvar, ncons
end


mutable struct WrappedListofconstraint{T}
    constraints::Array{Vector{T}}
    constants::Vector{T}
    islowerthan::BitVector
    isequalto::BitVector
    isgreaterthan::BitVector
    function Listofconstraint{T}()
        return new{T}([], [], falses(0), falses(0), falses(0))
    end
end
function numberofvarandconstraint(wlc::WrappedListofconstraint{T})::Tuple{Int, Int}
    nvar = sum(islowerthan) + sum(isgreaterthan)
    ncons = sum(islowerthan) + sum(isequalto) + sum(isgreaterthan)
    return nvar, ncons
end
mutable struct WrappedObjective{T}
    c::Vector{T}
    constant::T
    isminmax::Bool
    ismin::Bool
    function WrappedObjective{T}()
        return new{T}([], zero(T), true, true)
    end
end
function getAbc(wlv::WrappedListofvar{T}, wlc::WrappedListofconstraint{T}, wo::WrappedObjective{T}) where T
    nvar1, ncons1 = numberofvarandconstraint(wlv)
    nvar2, ncons2 = numberofvarandconstraint(wlc)
    nvar = nvar1 + nvar2
    ncons = ncons1 + ncons2
    A = zeros(T, ncons, nvar)
    b = zeros(T, ncons)
    c = zeros(T, nvar)

    currentvar = 1
    currentcons = 1
    for i in 