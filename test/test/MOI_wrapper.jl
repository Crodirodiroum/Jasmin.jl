# ============================ /test/MOI_wrapper.jl ============================
module TestWrap

import Jasmin
using MathOptInterface
using Test

#println("package imported and used successfuly")
const MOI = MathOptInterface

const OPTIMIZER = MOI.instantiate(
    MOI.OptimizerWithAttributes(Jasmin.Optimizer, MOI.Silent() => true),
)

const BRIDGED = MOI.instantiate(
    MOI.OptimizerWithAttributes(Jasmin.Optimizer, MOI.Silent() => true),
    with_bridge_type = Float64,
)

# See the docstring of MOI.Test.Config for other arguments.
const CONFIG = MOI.Test.Config(
    # Modify tolerances as necessary.
    atol = 1e-6,
    rtol = 1e-6,
    # Use MOI.LOCALLY_SOLVED for local solvers.
    optimal_status = MOI.OPTIMAL,
    # Pass attributes or MOI functions to `exclude` to skip tests that
    # rely on this functionality.
    exclude = Any[MOI.VariableName, MOI.delete],
)

#println("const defined")
"""
    runtests()

This function runs all functions in the this Module starting with `test_`.
"""
function runtests()
    #io = open("test_ran.txt", "w")
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            #write(io, "$(name)")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    #close(io)
end

"""
    test_runtests()

This function runs all the tests in MathOptInterface.Test.

Pass arguments to `exclude` to skip tests for functionality that is not
implemented or that your solver doesn't support.
"""
function test_runtests()
    MOI.Test.runtests(
        BRIDGED,
        CONFIG,
        exclude = String[ 
            "test_attribute_NumberOfThreads",#0/0
            "test_quadratic_",#0/0
            "test_conic",#129/217
            "test_add_constrained",#6/6
            "test_attribute", #9/9
            "test_basic", #289/289
            "test_constraint",#0/0
            "test_cpsat", #0/0
            #"test_linear", #189/283
            "test_model",#143/145 #no exception thrown
            "test_modification",# 107/164
            "test_nonlinear",#23/23
            "test_objective",#35/47
            "test_solve",#24/38
            "test_variable", #59/70 
            
        ],
        # This argument is useful to prevent tests from failing on future
        # releases of MOI that add new tests. Don't let this number get too far
        # behind the current MOI release though! You should periodically check
        # for new tests in order to fix bugs and implement new features.
        exclude_tests_after = v"0.10.5",
    )
    return
end


"""
    test_SolverName()

You can also write new tests for solver-specific functionality. Write each new
test as a function with a name beginning with `test_`.
"""
#=
function test_SolverName()
    @test MOI.get(Jasmin.Optimizer(), MOI.SolverName()) == "Jasmin"
    return
end
=#
end # module TestWrap

# This line at tne end of the file runs all the tests!
TestWrap.runtests()