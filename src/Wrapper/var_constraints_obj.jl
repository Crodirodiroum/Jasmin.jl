
#we only support greater than constraints for variables
function MOI.supports_constraint(opt::Optimizer, 
        ::Type{MOI.VariableIndex}, 
        ::Type{T}) where {T <: MOI.AbstractScalarSet}
    #println("constraint of type $T on variable index is not supported by the solver")
    return false
end
function MOI.supports_constraint(opt::Optimizer, 
        ::Type{MOI.VariableIndex}, 
        ::Type{MOI.GreaterThan{T}}) where T
    return true
end

#we only support affine function constraints that are equal to something as constriants.
function MOI.supports_constraint(
        opt::Optimizer,
        f::Type{F},
        s::Type{S},)::Bool where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    return false
end
function MOI.supports_constraint(
        opt::Optimizer,
        ::Type{MOI.ScalarAffineFunction{T}},
        ::Type{MOI.EqualTo{T}})::Bool where T
    return true
end


function MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {T, F <: MOI.ScalarAffineFunction{T}}
    return true
end
function MOI.supports(::Optimizer, ::MOI.ObjectiveSense)
    return true
end


function MOI.get(model::Optimizer,
        attr::MOI.ConstraintPrimal,
        c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},)
    #println("------------\n"^64)
    return 0.0#getvarvalue(c, model.map, model.vars, model.wrap.xstar)
end
function MOI.get(model::Optimizer,
        attr::MOI.ConstraintPrimal,
        c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},<:Any},) where T
    return zero(T)
end

function MOI.get(model::Optimizer, cd::MOI.ConstraintDual, ci::MOI.ConstraintIndex{<:Any, <:Any}) where T
    return zero(Float64)
end


function MOI.get(model::Optimizer, ::MOI.ListOfConstraintIndices{T, S}) where {T, S}
    return MOI.get(model.inner, MOI.ListOfConstraintIndices{T, S}())
end

function MOI.get(model::Optimizer, ::MOI.ConstraintFunction, ci::MOI.ConstraintIndex{T, S}) where {T, S}
    return MOI.get(model.inner, MOI.ConstraintFunction(), ci)
end