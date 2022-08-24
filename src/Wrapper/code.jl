mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    inner
    wrap::LinearProblem{T}
    isempty::Bool
    name::String
    silent::Bool
    attributs::Dict
    objsense::MOI.OptimizationSense
    vars::Array{WrappedVar{T}, 1}
    scalar_residus::T #something to add to the objective function
    vstar::T
    xstar::Vector{T}
    status::MOI.TerminationStatusCode
    timelimit::Float64
    solve_time::Float64
    map::MOI.IndexMap
    function Optimizer(::Type{T} = Float64) where T
        opt = new{T}()
        opt.isempty = true
        opt.name = ""
        opt.silent = true
        #opt.inner = MOI.FileFormats.NL.Model()
        opt.attributs = Dict()
        opt.timelimit = 5.0
        opt.vstar = zero(T)
        opt.xstar = zeros(T,0)
        opt.status = MOI.OPTIMIZE_NOT_CALLED
        opt.map = MOI.IndexMap()
        return opt
    end
end
#incremental disable
function MOI.supports_incremental_interface(::Optimizer)
    return false
end
#basic functions
function MOI.empty!(opt::Optimizer{T}) where T
    opt.isempty = true
    opt.wrap = LinearProblem{T}(zeros(T, 0, 0), zeros(T, 0), zeros(T, 0))
    opt.status = MOI.OPTIMIZE_NOT_CALLED
    return nothing
end
function MOI.is_empty(opt::Optimizer)
    return opt.isempty
end
function MOI.optimize!(opt::Optimizer{T}) where T
    tstart = time()
    inner = opt.wrap
    status = solve!(inner)
    #println("solve jasmin linearproblem just called")
    #println(inner)
    inner.vstar += opt.scalar_residus
    if opt.objsense == MOI.MAX_SENSE
        opt.vstar *= -1
    end
    opt.vstar = vstar(inner)
    opt.xstar = xstar(inner)
    opt.status = status
    opt.solve_time = time() - tstart
    return opt.status
end
function Base.show(io::IO, model::Optimizer)
    return print(io, "NewSolver $(typeof(model))")
end


#get
function MOI.get(optimizer::Optimizer, attr::MOI.SolverName)
    return "Jasmin"
end
function MOI.get(optimizer::Optimizer, attr::MOI.SolverVersion)
    #this is fixed by hand but should go check in the manifest for the package version of 
    return "v0.2.0"
end
function MOI.get(optimizer::Optimizer, attr::MOI.RawSolver)
    return optimizer.wrap
end
function MOI.get(optimizer::Optimizer, attr::MOI.Name)
    return optimizer.name
end
function MOI.get(optimizer::Optimizer, attr::MOI.Silent)
    return optimizer.Silent
end 
function MOI.get(optimizer::Optimizer, attr::MOI.TimeLimitSec)
    return optimizer.timelimit
end
function MOI.get(optimizer::Optimizer, attr::MOI.NumberOfThreads)
    #todo to multithreading stuff
    return 1
end

#set
function MOI.set(optimizer::Optimizer, attr::MOI.Name, value::String)
    optimizer.name = value
end
function MOI.set(optimizer::Optimizer, attr::MOI.Silent, value::Bool)
    optimizer.silent = value
end
function MOI.set(optimizer::Optimizer, attr::MOI.TimeLimitSec, t::Float64)
    optimizer.timelimit = t
end
function MOI.set(optimizer::Optimizer, attr::MOI.NumberOfThreads, n::Int)
    #todo
end

#supports
function MOI.supports(optimizer::Optimizer, attr::MOI.Name, value::String)
    return true
end
function MOI.supports(optimizer::Optimizer, attr::MOI.Silent, value::Bool)
    return true
end
function MOI.supports(optimizer::Optimizer, attr::MOI.TimeLimitSec)
    return true
end

function MOI.supports(optimizer::Optimizer, attr::MOI.NumberOfThreads)
    return false
end

function MOI.supports(model::Optimizer, ::MathOptInterface.RawStatusString)
    return true
end
function MOI.get(model::Optimizer, attr::MOI.RawStatusString)
    return "this should tell the user why the solver stop."
end



MOI.supports(::Optimizer{T}, ::MOI.RawOptimizerAttribute) where T = true

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    key = Symbol(param.name)
    model.attributs[key] = value
    println("raw attributs $key set to value")
    return
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    key = Symbol(param.name)
    println("raw attributs $key requested")
    return model.attributs[key]
end
#copy_to
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)::MOI.IndexMap
    dest.inner = src
    #println("cal to copy_to")
    T = Float64
    current = 1
    map = MOI.Utilities.IndexMap()
    vars = createWrappedVar!(src, map, current, T)
    dest.vars = vars
    A, b = createwrappedconstraint(src, map, current, vars, T)
    c = zeros(T, size(A, 2))
    dest.objsense = MOI.get(src, MOI.ObjectiveSense())
    if dest.objsense != MOI.FEASIBILITY_SENSE
        f = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
        row = createrow(f.terms, map, vars)
        c[1:length(row)] = row
        constant = f.constant
        if dest.objsense == MOI.MAX_SENSE
            c[:] *= -1
            constant *= -1
        end
        dest.scalar_residus = constant
    end
    jlp = LinearProblem{T}(A, b, c)
    jlp.timelimit = dest.timelimit
    dest.wrap = jlp
    #println(jlp)
    dest.map = map
    return map
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    return model.status
end

function MOI.get(model::Optimizer{T}, ::MOI.VariablePrimal, vi::MOI.VariableIndex)::T where T
    x = model.xstar
    #println("get value called")
    #getvarvalue(vi, map::MOI.Utilities.IndexMap, vars::Array{WrappedVar{T}}, x::Vector{T})
    getvarvalue(vi, model.map, model.vars, x)
end

function MOI.get(model::Optimizer, attr::MOI.SolveTimeSec)
    return model.solve_time
end


function MOI.get(model::Optimizer, ::MOI.PrimalStatus)
    status = model.status
    if status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif status == MOI.INFEASIBLE
        return MOI.NO_SOLUTION
    elseif status == MOI.DUAL_INFEASIBLE
        return MOI.NO_SOLUTION
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

MOI.supports(model::Optimizer, ::MOI.ObjectiveValue) = true
function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    return model.vstar
end

#return a MathOptInterface.ResultStatusCode
function MOI.get(::Optimizer, ::MOI.DualStatus)
    return MOI.NO_SOLUTION
end
function MOI.get(model::Optimizer, ::MOI.ResultCount)
    if model.status == MOI.OPTIMAL
        return 1
    else
        return 0
    end
end
function MOI.get(model::Optimizer{T}, ::MOI.DualObjectiveValue) where T
    model.status == MOI.OPTIMAL && return model.vstar
    model.objsense == MOI.MIN_SENSE && return Inf
    model.objsense == MOI.MAX_SENSE && return -Inf
    return zero(T)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
        c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any})
        #constraint on a single variable, return the value of the variable
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end




function MOI.get(model::Optimizer, ::MOI.VariableBasisStatus, ci::MOI.VariableIndex)::MOI.BasisStatusCode
    return MOI.NONBASIC_AT_LOWER
end
function MOI.get(model::Optimizer, ::MOI.ConstraintBasisStatus, ci::MOI.ConstraintIndex)::MOI.BasisStatusCode
    return MOI.BASIC
end
