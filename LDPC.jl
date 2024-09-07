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

NREALS::Int = 1_00
MAX::Int = 30

#################################### CODING ####################################

### Moreira's example (p.284)

# include("Moreira.jl")

### PEG

N = 512
M = 256
D = rand([2,3],N)

@time H, girth = PEG!(D,M)

println("girth = ", girth)

@time _, _, G = GF2_nullspace(H)

K = size(G,2)

println("K = ", K)

σ_test = 0.8

n_test = σ_test*randn(N)

Message = rand([0,1],K)

### codeword

C = G*Message .% 2

##################################### BPKS #####################################

U = 2*C .- 1


SNR = collect(3:-0.3:0.3)
Sigma = 1 ./ SNR

############################# AUXILIARY CONSTANTS ##############################

indices_M  = findindices_M(H)
indices_N  = findindices_N(H)

##################################### TEST #####################################
t_test = U + n_test
# RR, LRR, QQ, LQQ = test_SPA(indices_M,indices_N,Phi,t_test,σ_test,"TNH")

############################## JULIA COMPILATION ###############################

# performance_estimation(C,U,Sigma,H,indices_N,indices_M,Phi,"TNH";nreals=1)
performance_estimation(C,U,Sigma,H,indices_N,indices_M,Phi,"ALT";nreals=1)
performance_estimation(C,U,Sigma,H,indices_N,indices_M,Phi,"TAB";nreals=1)
performance_estimation(C,U,Sigma,H,indices_N,indices_M,Phi,"APP";nreals=1)
                             
########################### PERFORMANCE SIMULATION ############################

# @time FER_tnh, BER_tnh, Iters_tnh = performance_estimation(C,U,Sigma,H,indices_N,
#                                                             indices_M,Phi,"TNH")
@time FER_alt, BER_alt, Iters_alt = performance_estimation(C,U,Sigma,H,indices_N,
                                                            indices_M,Phi,"ALT")
@time FER_tab, BER_tab, Iters_tab = performance_estimation(C,U,Sigma,H,indices_N,
                                                            indices_M,Phi,"TAB")
@time FER_app, BER_app, Iters_app = performance_estimation(C,U,Sigma,H,indices_N,
                                                          indices_M, Phi, "APP")

################################### PLOTTING ###################################
plotlyjs()
lim = log10(1/NREALS)
# yaxis = [FER_tnh, FER_alt ,FER_tab, FER_app]
# labels = permutedims(["SPA 1", "SPA 2", "Lookup Table SPA", "Min Sum"])
# display(plot(SNR, yaxis, label=labels, linewidth=2, title="FER", ylims=(lim,0)))

yaxis = [FER_alt ,FER_tab, FER_app]
labels = permutedims(["SPA", "Lookup Table SPA", "Min Sum"])
display(plot(SNR, yaxis, label=labels, linewidth=2, title="FER", ylims=(lim,0)))

plot(1:MAX, BER_alt, linewidth=2, title="BER SPA", ylims=(lim-1,0))
plot!(girth*[1, 1]/2, [lim-1, 0], linewidth=2, linestyle=:dot, linecolor=:black)

