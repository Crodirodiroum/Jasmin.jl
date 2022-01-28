include("problems.jl")


function testProblem(::Type{so}, ::Type{l}, p::NamedTuple, ::Type{T} = Float64)::Bool where {so <: AbstractSolver, l <: AbstractLP, T}
    (sol, vsol,  A, b, c, isoptimal, isbounded, isfeasible) = (p.sol, p.vsol,  p.A, p.b, p.c, p.isoptimal, p.isbounded, p.isfeasible)
    linearform = l{T}(A, b, c)
    changeSolver!(so)
    solve!(linearform)
    goodstatus = false
    if isoptimal
        goodstatus = (linearform.status == Optimal)
    elseif !isbounded
        goodstatus = (linearform.status == Unbounded)
        return true
    elseif !isfeasible
        goodstatus = (linearform.status == Infeasible)
        return true
    end

    goodsol = isapprox(minimum([norm(xstar(linearform) - s) for s in sol]), 0, atol = 1e-6)
    goodsolvalue = isapprox(vstar(linearform), vsol, atol = 1e-6)
    #@show goodsol
    #@show goodsolvalue
    #@show goodstatus
    #@show linearform.status
    return goodsol && goodsolvalue && goodstatus
end

solverlist = [BaseEnumeration, StandardSimplexe]
representationlist = [LinearProblem]


for so in solverlist
    for rep in representationlist
        for p in ProblemTest.PL
            @test testProblem(so, rep, p)
        end
    end
end
