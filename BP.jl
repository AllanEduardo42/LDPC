using LinearAlgebra
using Statistics

include("horizontal_update.jl")
include("simple_horizontal_update.jl")
include("vertical_update.jl")
include("MAP_estimator.jl")
include("auxiliary_functions.jl")
include("log_simple_horizontal_update.jl")
include("log_vertical_update.jl")
include("log_MAP_estimator.jl")
include("llr_simple_horizontal_update.jl")
include("llr_vertical_update.jl")
include("llr_MAP_estimator.jl")
include("lookUPtable.jl")

SIZE = 512
SCOPE = SIZE/10
SCOPE_LN2 = round(Int,SCOPE*log(2))
lookUPtable()

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

t = 2*c .- 1

σ = 0.8

n = σ*randn(N)

r = t .+ n

r = [1.3129, 2.6584, 0.7413, 2.1745, 0.5981, −0.8323, −0.3962, −1.7586,
1.4905, 0.4084, −0.9290, 1.0765]

f = zeros(N,2)
k = 1/(sqrt(2π)*σ)
s = 2σ^2
for i in eachindex(r)
    f[i,1] = k*exp(-(r[i]+1)^2/s)
    f[i,2] = k*exp(-(r[i]-1)^2/s)
end

normalize(f,N)

LF = abs.(log.(f))

indices_M  = findindices_M(H,N)
indices_N  = findindices_N(H,M)

Q = init_q(M,N,f,indices_M)
LQ = log_init_q(M, N, LF, indices_M)
LLQ = llr_init_q(M, N, LF, indices_M)

for i=1:3

    # global R = horizontal_update(M, N, Q, indices_N)

    @time global R = simple_horizontal_update(M, N, Q, indices_N)

    @time global LR = log_simple_horizontal_update(M,N,LQ,indices_N, indices_M)

    @time global LLR = llr_simple_horizontal_update(M,N,LLQ,indices_N, indices_M)

    global Q = vertical_update(N,R,Q,f, indices_M)

    global LQ = log_vertical_update(N, LR, LQ, LF, indices_M)

    global LLQ = llr_vertical_update(N, LLR, LLQ, LF, indices_M)

end

d = MAP_estimator(N,R,f,indices_M)

Ld = log_MAP_estimator(N, LR, LF, indices_M)

LLd = llr_MAP_estimator(N, LLR, LF, indices_M)

println(rem.(H*d,2))
println(rem.(H*Ld,2))
println(rem.(H*LLd,2))