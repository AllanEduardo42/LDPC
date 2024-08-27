using LinearAlgebra
using Random
using Statistics
using Plots

plotlyjs()

# include("horizontal_update.jl")
include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")
include("auxiliary_functions.jl")
include("llr_simple_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")
include("performance_estimation.jl")

NREALS::Int = 1000
MAX::Int = 10

H = Bool.([0 1 0 1 0 1 1 1 0 0 0 1;
     1 0 1 1 0 0 0 0 1 0 0 0;
     0 1 0 0 1 0 1 0 0 0 0 1;
     1 0 0 1 0 0 0 0 0 1 1 0;
     0 0 1 0 1 1 0 0 0 1 0 0;
     1 0 1 0 0 0 1 1 0 0 1 0;
     0 1 0 0 0 1 0 1 1 1 0 0;
     0 0 0 0 1 0 0 0 1 0 1 1])

MM,NN = size(H)
K = NN - MM

A = H[:,1:MM]
B = H[:,MM+1:NN]

P = Bool.(abs.(rem.(A\B,2)))

G = [P; I(K)]

Message = Bool.([1, 0, 0, 0])

C = rem.(G*Message,2)

R = 2*C .- 1

Sigma = [0.8]

indices_M  = findindices_M(H,NN)
indices_N  = findindices_N(H,MM)

performance_estimation(R, Sigma, MM, NN, indices_N, indices_M, H; nreals = 1)

# @time FER, Iters, LLQ, LR = performance_estimation(R, Sigma, MM, NN, indices_N, indices_M, H)

plot(log10.(FER[end:-1:1]))