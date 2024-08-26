using LinearAlgebra
using Statistics

# include("horizontal_update.jl")
include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")
include("auxiliary_functions.jl")
include("llr_simple_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")

H = Bool.([0 1 0 1 0 1 1 1 0 0 0 1;
     1 0 1 1 0 0 0 0 1 0 0 0;
     0 1 0 0 1 0 1 0 0 0 0 1;
     1 0 0 1 0 0 0 0 0 1 1 0;
     0 0 1 0 1 1 0 0 0 1 0 0;
     1 0 1 0 0 0 1 1 0 0 1 0;
     0 1 0 0 0 1 0 1 1 1 0 0;
     0 0 0 0 1 0 0 0 1 0 1 1])

M,N = size(H)
K = N - M

# Ht = [I_M P]

A = H[:,1:M]
B = H[:,M+1:N]

P = Bool.(abs.(rem.(A\B,2)))

G = [P; I(K)]

m = Bool.([1, 0, 0, 0])

c = Bool.(rem.(G*m,2))

σ = 0.8

# n = σ*randn(N)

# t = (2*c .- 1 ).+ n

t = [1.3129, 2.6584, 0.7413, 2.1745, 0.5981, −0.8323, −0.3962, −1.7586,
1.4905, 0.4084, −0.9290, 1.0765]

F = zeros(N,2)
k = 1/(sqrt(2π)*σ)
s = 2σ^2
for i in eachindex(t)
    F[i,1] = k*exp(-(t[i]+1)^2/s)
    F[i,2] = k*exp(-(t[i]-1)^2/s)
end

normalize(F,N)

ΔLF = log.(F[:,1]) - log.(F[:,2])

indices_M  = findindices_M(H,N)
indices_N  = findindices_N(H,M)

Q = init_q(M,N,F,indices_M)
LQ = llr_init_q(M, N, ΔLF, indices_M)

S = -1
i = 0
MAX = 10

while S != 0 && i < MAX

    global i += 1

    println("Iteration #", i)

    println("Conventional SPA:")

    @time global R = simple_horizontal_update(M, N, Q, indices_N)
    @time global Q, D = vertical_update_and_MAP(N,R,Q,F, indices_M)

    syndrome = rem.(H*D,2)

    println("MAP estimate syndrome: ", syndrome)

    println("LLR SPA:")

    @time global LR = llr_simple_horizontal_update(M,N,LQ,indices_N)
    @time global LQ, D_llr = llr_vertical_update_and_MAP(N, LR, LQ, ΔLF, indices_M)

    syndrome_llr = rem.(H*D_llr,2)    

    println("LLR MAP estimate syndrome: ", syndrome_llr)
    
    global S = sum(syndrome_llr)

end