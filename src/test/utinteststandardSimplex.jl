const PT = ProblemTest

changeSolver!(StandardSimplexe)
lp1Float64 = LinearProblem{Float64}(PT.P1.A, PT.P1.b, PT.P1.c)
solve!(lp1Float64, verbose = false)
@assert minimum([norm(xstar(lp1Float64) - sol) for sol in PT.P1.sol]) <= 1e-6 "wrong solution"
@assert vstar(lp1Float64) == -1.0 "wrong solution"
@assert typeof(vstar(lp1Float64)) == Float64 "type lost"

lp1RationalInt64 = LinearProblem{Rational{Int64}}(PT.P1.A, PT.P1.b, PT.P1.c)
solve!(lp1RationalInt64)
#test trivial problem but with Rational{Int}
@assert minimum([norm(xstar(lp1RationalInt64) - sol) for sol in PT.P1.sol]) <= 1e-6 "wrong solution"
@assert vstar(lp1RationalInt64) == PT.P1.vsol "wrong solution"
@assert typeof(vstar(lp1RationalInt64)) == Rational{Int64} "type lost"

lp2 = LinearProblem{Float64}(PT.P2.A, PT.P2.b, PT.P2.c)
solve!(lp2, verbose = false)
#test trivial problem with duplicated constraint
@assert minimum([norm(xstar(lp2) - sol) for sol in PT.P2.sol]) <= 1e-6 "wrong solution"
@assert vstar(lp2) == PT.P2.vsol "wrong solution"
@assert typeof(vstar(lp2)) == Float64 "type lost"

for pkmk in PT.PKM
    A, b, c, vs, xs = (pkmk.A, pkmk.b, pkmk.c, pkmk.vsol,pkmk.sol)
    lp = LinearProblem{Float64}(A, b, c)
    solve!(lp, verbose = false)
    @assert norm(xstar(lp) - xs) <= 1e-9 "wrong solution for Klee Minty no $k"
    @assert norm(vstar(lp) - vs) <= 1e-9 "wrong solution for Klee Minty no $k"
end