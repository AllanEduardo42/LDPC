################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Compare residue decaying for RBP
################################ JULIA PACKAGES ################################

using LinearAlgebra
using Random
using Statistics
using Plots
using SparseArrays
using CSV, DataFrames

################################ INCLUDED FILES ################################

include("performance_simulation.jl")
include("PEG.jl")
include("GF2_functions.jl")

############################# SIMULATION CONSTANTS #############################

const INF = typemax(Int64)

# Seeds
SEED_NOISE::Int64 = 1428
SEED_GRAPH::Int64 = 5714
SEED_MESSA::Int64 = 9999

###################### NUMBER OF TRIALS AND MULTITHREADING #####################

TRIALS::Int = 32768
NTHREADS::Int = min(32,TRIALS)

######################## MAXIMUM NUMBER OF BP ITERATIONS #######################

MAXRBP::Int = 5
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are 
# not printed)

################################ CONTROL FLAGS #################################

SAVEDATA::Bool = false
PRINTTEST::Bool = false

################################# RBP CONSTANTS ################################

# residue decays
RD = collect(0.1:0.1:1)
DECAYRBP::Float64 = RD[end]

##################################### SNR ######################################
SNRTEST = [3]
SNR = collect(1:1:4)

############################# PARITY-CHECK MATRIX #############################

# Matrix dimensions
N::Int64 = 512
M::Int64 = 256

# Vector of the variable node degrees
D = rand(Xoshiro(SEED_GRAPH),[2,3,4],N)

# Generate Parity-Check Matrix by the PEG algorithm
H, girth = PEG(D,M)

# Find the generator matrix
G = gf2_nullspace(H)

############################# MESSAGE AND CODEWORD #############################

Message = rand(Xoshiro(SEED_MESSA),Bool,N-M)

Codeword = gf2_mat_mult(Matrix(G), Message)

######################### PRINT INFORMATION ON SCREEN ##########################
println()
print("############################### LDPC parameters #######################")
println("#########")
println()
println("Parity Check Matrix: $M x $N")
println()
display(sparse(H))
println()
println("Graph girth = ", girth)
println()
println("Message:")
for i in eachindex(Message)
    print(Int(Message[i]))
    if i%80 == 0
        println()
    end
end
println()
println()
println("Codeword:")
for i in eachindex(Codeword)
    print(Int(Codeword[i]))
    if i%80 == 0
        println()
    end
end
println()
println()

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

rgn_noise_seeds = zeros(Int,NTHREADS)
for i in eachindex(rgn_noise_seeds)
    rgn_noise_seeds[i] = SEED_NOISE + i - 1
end

#################### JULIA COMPILATION (FOR SPEED) AND TEST ####################
Lr , Lq = performance_simulation(Codeword,SNRTEST,H,"RBP",1,MAXRBP,STOP,
                                rgn_noise_seeds;printtest=PRINTTEST)
                           
############################ PERFORMANCE SIMULATION ############################
if TRIALS > 1
    fer_labels = Vector{String}()
    fermax = Vector{Vector{<:AbstractFloat}}()
    FER = zeros(MAXRBP,length(SNR),length(RD))
    BER = zeros(MAXRBP,length(SNR),length(RD))
    for i in eachindex(RD)
            global DECAYRBP = RD[i]
            @time FER[:,:,i], BER[:,:,i] = performance_simulation(Codeword,
                                            SNR,H,"RBP",TRIALS,MAXRBP,STOP,
                                            rgn_noise_seeds)  
        push!(fer_labels,"R. decay = $(RD[i])")
        push!(fermax,FER[MAXRBP,:,i])      
    end
end
################################### PLOTTING ###################################
if TRIALS > 1
    plotlyjs()
    lim = log10(1/TRIALS)

    # FER x SNR
    if length(SNR) > 1
        fer_labels = permutedims(fer_labels)    
        p = plot(
            SNR,fermax,
            xlabel="SNR (dB)",
            label=fer_labels,
            lw=2,
            title="FER (Graph girth = $girth)",
            ylims=(lim,0)
        )
        display(p)
        SAVEDATA ? savefig(p, "FER.png") : nothing
    end
end