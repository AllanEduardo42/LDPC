################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# LDPC coding-decoding performance simulation

################################ JULIA PACKAGES ################################

using LinearAlgebra
using Random
using Statistics
using Plots

################################ INCLUDED FILES ################################

# include("horizontal_update.jl")
include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")
include("auxiliary_functions.jl")
include("llr_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")
include("performance_estimation.jl")
include("test_SPA.jl")
include("lookupTable.jl")

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1428

SIZE::Int64 = 1024
RANGE::Int64 = 10

SIZE_per_RANGE::Float64 = SIZE/RANGE

Phi = lookupTable()

NREALS::Int = 100000
MAX::Int = 10

#################################### CODING ####################################

H = [0 1 0 1 0 1 1 1 0 0 0 1;
     1 0 1 1 0 0 0 0 1 0 0 0;
     0 1 0 0 1 0 1 0 0 0 0 1;
     1 0 0 1 0 0 0 0 0 1 1 0;
     0 0 1 0 1 1 0 0 0 1 0 0;
     1 0 1 0 0 0 1 1 0 0 1 0;
     0 1 0 0 0 1 0 1 1 1 0 0;
     0 0 0 0 1 0 0 0 1 0 1 1]

MM,NN = size(H)
K = NN - MM

A = H[:,1:MM]
B = H[:,MM+1:NN]

P = abs.(Int64.(A\B .% 2))

G = [P; I(K)]

Message = [1, 0, 0, 0]

C = G*Message .% 2

#################################### BPKS ####################################

U = 2*C .- 1

Sigma = collect(0.1:0.1:1)

############################# AUXILIARY CONSTANTS #############################

indices_M  = findindices_M(H,NN)
indices_N  = findindices_N(H,MM)

#################################### TEST #####################################

RR, LRR, QQ, LQQ = test_SPA(indices_M, indices_N, H, MM, NN, Phi)

############################# JULIA COMPILATION ##############################

___ = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, H;
                             nreals = 1)
___ = performance_estimation_table(C, U, Sigma, MM, NN, indices_N, indices_M, H, Phi;
                              nreals = 1)

########################### PERFORMANCE SIMULATION ############################

@time FER, Iters, ___, ___ = performance_estimation(C, U, Sigma, MM, NN,
                                                    indices_N, indices_M, H)
@time FERt, Iterst, ___, ___ = performance_estimation_table(C, U, Sigma, MM, NN,
                                                    indices_N, indices_M, H, Phi)

plotlyjs()

################################### PLOTTING ###################################

plot(log10.(FERt[end:-1:1]))
plot!(log10.(FER[end:-1:1]))
plot!([1,10],log10.(1/NREALS*[1, 1]))