using LinearAlgebra
using SparseArrays
using Random

include("GF2_functions.jl")

function encoding_1(
    trials::Int,
    H::Matrix{Bool},
    L::Matrix{Bool},
    U::Matrix{Bool}
)

    rgn = Xoshiro(1234)

    M,N = size(H)
    K = N - M

    r = zeros(Bool,M)
    cword = zeros(Bool,N)
    c = zeros(Bool,M)
    H1 = H[:,1:K]

    for i in 1:trials
        rand!(rgn,msg,Bool)
        c = H1*msg
        r = gf2_solve_LU(L,U,c)
        cword[1:K] .= msg
        cword[K+1:end] .= r
    end

end

function encoding_2(
    trials::Int,
    H::Matrix{Bool},
    L::Matrix{Bool},
    U::Matrix{Bool}
)

    M,N = size(H)
    K = N - M

    P = zeros(Bool,M,K)
    H1 = H[:,1:K]
    r = zeros(Bool,M)
    cword = zeros(Bool,N)

    for k in axes(P,2)
        P[:,k] = gf2_solve_LU(L,U,H1[:,k])
    end

    for i in 1:trials
        rand!(rgn,msg,Bool)
        r = P*msg
        cword[1:K] .= msg
        cword[K+1:end] .= r
    end

end

trials = 10000

@benchmark encoding_1($trials,$H,$L,$U)

@benchmark encoding_2($trials,$H,$L,$U)

