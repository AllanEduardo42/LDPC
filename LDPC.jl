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
include("SPA.jl")
include("PEG.jl")
include("GF2_nullspace.jl")

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1427

SIZE::Int64 = 1024
RANGE::Int64 = 20

SIZE_per_RANGE::Float64 = SIZE/RANGE

Phi = lookupTable()

NREALS::Int = 2_000_000
MAX::Int = 3

#################################### CODING ####################################

### PEG

NN = 60
MM = 30
# d = rand([2,3],NN)
# sort!(d)
# d = 2*ones(Int,NN)
d = 3*ones(Int,NN)

H, girth = PEG(MM,NN,d)

println("grith = ", girth)

_, _, G = GF2_nullspace(H)

K = size(G,2)

println("K = ", K)

### Moreira's example using Mckay method:

# H = [0 1 0 1 0 1 1 1 0 0 0 1;
#      1 0 1 1 0 0 0 0 1 0 0 0;
#      0 1 0 0 1 0 1 0 0 0 0 1;
#      1 0 0 1 0 0 0 0 0 1 1 0;
#      0 0 1 0 1 1 0 0 0 1 0 0;
#      1 0 1 0 0 0 1 1 0 0 1 0;
#      0 1 0 0 0 1 0 1 1 1 0 0;
#      0 0 0 0 1 0 0 0 1 0 1 1]

# MM,NN = size(H)
# K = NN - MM

# A = H[:,1:MM]
# B = H[:,MM+1:NN]

# P = abs.(Int64.(A\B .% 2))

# G = [P; I(K)]

### Message

Message = ones(Int,K)

C = G*Message .% 2

##################################### BPKS #####################################

U = 2*C .- 1


SNR = collect(3:-0.3:0.3)
Sigma = 1 ./ SNR

############################# AUXILIARY CONSTANTS ##############################

indices_M  = findindices_M(H,NN)
indices_N  = findindices_N(H,MM)

##################################### TEST #####################################

RR, LRR, QQ, LQQ = test_SPA(indices_M,indices_N,H,MM,NN,Phi,"APP")

############################## JULIA COMPILATION ###############################


performance_estimation(C,U,Sigma,MM,NN,indices_N,indices_M,Phi,"TNH";nreals=1)
performance_estimation(C,U,Sigma,MM,NN,indices_N,indices_M,Phi,"ALT";nreals=1)
performance_estimation(C,U,Sigma,MM,NN,indices_N,indices_M,Phi,"TAB";nreals=1)
performance_estimation(C,U,Sigma,MM,NN,indices_N,indices_M,Phi,"APP";nreals=1)
                             
########################### PERFORMANCE SIMULATION ############################

@time FER_tnh, Iters_tnh = performance_estimation(C,U,Sigma,MM,NN,indices_N,
                                                            indices_M,Phi,"TNH")
@time FER_alt, Iters_alt = performance_estimation(C,U,Sigma,MM,NN,indices_N,
                                                            indices_M,Phi,"ALT")
@time FER_tab, Iters_tab = performance_estimation(C,U,Sigma,MM,NN,indices_N,
                                                            indices_M,Phi,"TAB")
@time FER_app, Iters_app = performance_estimation(C,U,Sigma,MM,NN,indices_N,
                                                          indices_M, Phi, "APP")

################################### PLOTTING ###################################
plotlyjs()
lim = log10(1/NREALS)
yaxis = [FER_tnh ,FER_alt ,FER_tab, FER_app]
labels = permutedims(["Tanh", "Alt", "Table", "Approx"])
plot(SNR, yaxis, label=labels, linewidth=2, title="FER", ylims=(lim,0))
