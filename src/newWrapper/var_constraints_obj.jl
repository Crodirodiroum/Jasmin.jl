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


function MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: MOI.AbstractScalarFunction}
    return true
end
function MOI.supports(::Optimizer, ::MOI.ObjectiveSense)
    return true
end

