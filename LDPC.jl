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

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1428

SIZE::Int64 = 1024
RANGE::Int64 = 20

SIZE_per_RANGE::Float64 = SIZE/RANGE

Phi = lookupTable()

NREALS::Int = 1_000_00
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

RR, LRR, QQ, LQQ = test_SPA(indices_M, indices_N, H, MM, NN, Phi, "APPROX")

############################# JULIA COMPILATION ##############################


_ = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "TANH"; nreals = 1)
_ = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "ALT"; nreals = 1)
_ = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "TABLE"; nreals = 1)
_ = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "APPROX"; nreals = 1)
                             
########################### PERFORMANCE SIMULATION ############################

@time FER_tanh, Iters_tanh = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "TANH")
@time FER_alt, Iters_alt = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "ALT")
@time FER_table, Iters_table = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "TABLE")
@time FER_approx, Iters_approx = performance_estimation(C, U, Sigma, MM, NN, indices_N, indices_M, Phi, "APPROX")

plotlyjs()

################################### PLOTTING ###################################
yaxis = [log10.(FER_tanh), log10.(FER_alt), log10.(FER_table), log10.(FER_approx)]
labels = permutedims(["Tanh", "Alt", "Table", "Approx"])
p = plot(Sigma, yaxis, label = labels, linewidth = 2, title = "FER")
display(p)
plot!([minimum(Sigma), maximum(Sigma)],log10.(1/NREALS*[1, 1]),label="min")
# plot(Sigma,log10.(FER_tanh))
# plot!(Sigma,log10.(FER_alt))
# plot!(Sigma,log10.(FER_table))
# plot!(Sigma,log10.(FER_approx))
