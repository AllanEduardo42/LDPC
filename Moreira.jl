###############################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Test LDPC using Moreira's example

##################################### TEST #####################################

using LinearAlgebra

include("GF2_functions.jl")
include("auxiliary_functions.jl")
include("performance_estimation.jl")

PRINTING = true
SEED = 1428

MAX = 3
NREALS = 1

SIZE::Int64 = 1024
RANGE::Float64 = 20
SIZE_per_RANGE::Float64 = SIZE/RANGE

Phi = lookupTable()

mode = "MIN"

H = BitMatrix(
     [0 1 0 1 0 1 1 1 0 0 0 1;
      1 0 1 1 0 0 0 0 1 0 0 0;
      0 1 0 0 1 0 1 0 0 0 0 1;
      1 0 0 1 0 0 0 0 0 1 1 0;
      0 0 1 0 1 1 0 0 0 1 0 0;
      1 0 1 0 0 0 1 1 0 0 1 0;
      0 1 0 0 0 1 0 1 1 1 0 0;
      0 0 0 0 1 0 0 0 1 0 1 1]
    )

Nodes2checks  = find_nodes2checks(H)
Checks2nodes  = find_checks2nodes(H)

M,N = size(H)
K = N - M

### McKay Method

A = H[:,1:M]
B = H[:,M+1:N]

P = gf2_mat_mult(gf2_inverse(A),B)

G = [P; I(K)]

Message = Vector(Bool.([1, 0, 0, 0]))

C = gf2_mat_mult(Matrix(G), Message)

σ = [0.8]

t = [1.3129, 2.6584, 0.7413, 2.1745, 0.5981, −0.8323, −0.3962, −1.7586,
1.4905, 0.4084, −0.9290, 1.0765]

@time R, LR, Q, LQQ = 
    performance_estimation(
        C,
        σ,
        H,
        Checks2nodes,
        Nodes2checks,
        mode,
        1;
        t_test=t,
        printing=PRINTING
    )
;
if mode == "TAB"
    LR /= SIZE_per_RANGE
    LQQ /= SIZE_per_RANGE
end

logQ = log.(Q[:,:,1]./Q[:,:,2])
logR = log.(R[:,:,1]./R[:,:,2]);
