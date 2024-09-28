###############################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Test LDPC using Moreira's example

##################################### TEST #####################################

using LinearAlgebra

include("GF2_functions.jl")
include("auxiliary_functions.jl")
include("performance_simulation.jl")

PRINTING::Bool = true
SEED::Int = 1

MAX::Int = 3
NREALS::Int = 1

SIZE::Int = 1024
RANGE::Int = 20
SIZE_per_RANGE::Float64 = SIZE/RANGE

MODE::String = "MKAY"

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

Lr, Lq = 
    performance_simulation(C,σ,H,MODE,1,MAX;t_test=t,printing=PRINTING)
;
if MODE == "FTAB"
    Lr /= SIZE_per_RANGE
    Lq /= SIZE_per_RANGE
end

if MODE == "MKAY"
    logR = log.(Lr[:,:,1]./Lr[:,:,2])
end
