mutable struct WrappedVar{T}
    upper::T
    lower::T
    hasupperbound::Bool
    haslowerbound::Bool
    indexs::Tuple{Int, Int}
    vi::MOI.VariableIndex
    function WrappedVar{T}(vi::MOI.VariableIndex, current::Int) where T
        upper = T(Inf)
        lower = -T(Inf)
        hasupperbound = false
        haslowerbound = false
        indexs = (current, current + 1)
        vi = vi
        return new{T}(upper, lower, hasupperbound, haslowerbound, indexs, vi)
    end
end
#create wraped variable for all variables in the model, add them to the map and add all variable constraints to variables.
function createWrappedVar!(src::MOI.ModelLike, map::MOI.Utilities.IndexMap, current::Int, ::Type{T} = Float64) where T
    vars = MOI.get(src, MOI.ListOfVariableIndices())
    wvs = Array{WrappedVar{T}, 1}()
    secondcurrent = 1
    for v in vars
        wv = WrappedVar{T}(v, current)
        map[v] = MOI.VariableIndex(secondcurrent)
        secondcurrent+=1
        current += 2
        push!(wvs, wv)
    end
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.LessThan{T}}())
        varindex = MOI.get(src, MOI.ConstraintFunction(), ci)
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        wv = wvs[map[varindex].value]
        wv.upper = min(wv.upper, s.upper)
        wv.hasupperbound = true
    end
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.GreaterThan{T}}())
        varindex = MOI.get(src, MOI.ConstraintFunction(), ci)
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        wv = wvs[map[varindex].value]
        wv.lower = max(wv.lower, s.lower)
        wv.haslowerbound = true
    end
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.EqualTo{T}}())
        varindex = MOI.get(src, MOI.ConstraintFunction(), ci)
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        wv = wvs[map[varindex].value]
        wv.lower = max(wv.lower, s.lower)
        wv.upper = min(wv.upper, s.upper)
        wv.hasupperbound = true
        wv.haslowerbound = true
    end
    return wvs
end
function createrow(terms::Array{MathOptInterface.ScalarAffineTerm{T}, 1}, map::MOI.Utilities.IndexMap, vars::Array{WrappedVar{T}}) where T
    row = zeros(T, 2*length(vars))
    for term in terms
        coef, var = term.coefficient, term.variable
        var = vars[map[var].value]
        id1, id2 = var.indexs
        row[id1] = coef
        row[id2] = -coef
    end
    return row
end
function addrow(A::Matrix{T}) where T
    m, n = size(A)
    return [A; zeros(T, 1, n)]
end
function addcol(A::Matrix{T}) where T
    m, n = size(A)
    return [A zeros(T, m, 1)]
end
function createwrappedconstraint(src::MOI.ModelLike, map::MOI.Utilities.IndexMap, 
        current::Int, vars::Array{WrappedVar{T}}, ::Type{T} = Float64)::Tuple{Matrix{T}, Vector{T}} where T
    F = MOI.ScalarAffineFunction{T}
    A = zeros(T, 0, 2*length(vars))
    b = zeros(T, 0)
    nconstraint = 1
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F, MOI.EqualTo{T}}())
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        #println(prod("-" for _ in 1:80))
        #println("in $(MOI.ListOfConstraintIndices{F, MOI.EqualTo{T}}), $f")
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        row = createrow(f.terms, map, vars)
        coef = f.constant
        b = [b; set.value - coef]
        A = [A; row']
        map[ci] = MOI.ConstraintIndex{F, MOI.EqualTo{T}}(nconstraint)
        nconstraint += 1
    end
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F, MOI.LessThan{T}}())
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        #println(prod("-" for _ in 1:80))
        #println("in $(MOI.ListOfConstraintIndices{F, MOI.LessThan{T}}), $f")
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        row = createrow(f.terms, map, vars)
        coef = f.constant
        b = [b; set.upper - coef]
        m, n = size(A)
        A = addrow(A)
        A[end, 1:length(row)] = row
        A = addcol(A)
        current += 1
        A[end, end] = one(T)
        map[ci] = MOI.ConstraintIndex{F, MOI.LessThan{T}}(nconstraint)
        nconstraint += 1
    end
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F, MOI.GreaterThan{T}}())
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        #println(prod("-" for _ in 1:80))
        #println("in $(MOI.ListOfConstraintIndices{F, MOI.GreaterThan{T}}), $f")
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        row = createrow(f.terms, map, vars)
        coef = f.constant
        b = [b; set.upper - coef]
        m, n = size(A)
        A = addrow(A)
        A[end, 1:length(row)] = row
        A = addcol(A)
        current += 1
        A[end, end] = -one(T)
        
        map[ci] = MOI.ConstraintIndex{F, MOI.GreaterThan{T}}(nconstraint)
        nconstraint += 1
    end
    for var in vars
        if var.haslowerbound && var.hasupperbound && (var.lower == var.upper)
            A = addrow(A)
            id1, id2 = var.indexs
            A[end, id1] = one(T)
            A[end, id2] = -one(T)
            b = [b; var.lower]
        else
            if var.haslowerbound
                A = addrow(A)
                A = addcol(A)
                current += 1
                id1, id2 = var.indexs
                A[end, id1] = one(T)
                A[end, id2] = -one(T)
                A[end, end] = -1
                b = [b; var.lower]
            end
            if var.hasupperbound
                A = addrow(A)
                A = addcol(A)
                current += 1
                id1, id2 = var.indexs
                A[end, id1] = one(T)
                A[end, id2] = -one(T)
                A[end, end] = 1
                b = [b; var.upper]
            end
        end
    end
    return A, b
end

function getvarvalue(vi::MOI.VariableIndex, map::MOI.Utilities.IndexMap, vars::Array{WrappedVar{T}}, x::Vector{T})::T where T
    id1, id2 = vars[map[vi].value].indexs
    return x[id1] - x[id2]
end