module ProblemTest
using LinearAlgebra
P1 = (sol = [[1,0], [0, 1]], vsol = -1,  A = [1 1], b = [1], c = [-1, -1], isoptimal = true, isbounded = true, isfeasible = true)
P2 = (sol = [[1,0], [0, 1]], vsol = -1, A = [1 1; 1 1], b = [1, 1], c = [-1, -1], isoptimal = true, isbounded = true, isfeasible = true)
#=
function kmstar(n::Int, ::Type{T} = Int) where T
    xstar = zeros(T, 2*n)
    xstar[n] = 100^(n - 1)
    rest = [100^(k - 1) for k in 1:n]
    rest[end] = 0
    xstar[n+1:end] = rest
    return xstar
end
function km(n::Int, ::Type{T} = Int) where T
    c = [10^(n - i) for i in 1:n]
    b = [100^(i - 1) for i in 1:n]
    A = zeros(T, n, n)
    for i in 1:n
        A[i, 1:i] = 2*[10^(i - j) for j in 1:i]
    end
    A[:, :] -= I
    A = [A I]
    c = [c; zeros(T, n)]
    return A, b, -c, -100^(n - 1), kmstar(n)
end
PKM = [
    ((A, b, c, vstar, xstar) = km(i);
        (sol = [xstar],
        vsol = vstar,
            A = A,
            b = b,
            c = c, 
            isoptimal = true, 
            isbounded = true, 
            isfeasible = true
        )
    )
    for i in 1:10
]
=#
#Problem List
#PL = [P1, P2, PKM...]
PL = [P1, P2]#, PKM...]
end