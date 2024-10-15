################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# LDPC coding-decoding performance simulation using Believe Propagation (BP)

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
SEED_SAMPL::Int64 = 2857
SEED_MESSA::Int64 = 9999

###################### NUMBER OF TRIALS AND MULTITHREADING #####################

TRIALS::Int = 1024
NTHREADS::Int = min(32,TRIALS)

######################## MAXIMUM NUMBER OF BP ITERATIONS #######################

MAX::Int = 20
MAXRBP::Int = 5
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are 
# not printed)

################################ BP MODE FLAGS ################################

#Flooding
FLOO::Bool = false      
#LBP
_LBP::Bool = false      
#instantaneos-LBP
iLBP::Bool = false      
#RBP
_RBP::Bool = true      
#Random-RBP
RRBP::Bool = true      
#Local-RBP
LRBP::Bool = true      
#List-RBP
LIST::Bool = true      

# Aggregate mode string name and number of iterations
modes = [(FLOO,"Flooding",MAX),
         (_LBP,"LBP",MAX),
         (iLBP,"iLBP",MAX),
         (_RBP,"RBP",MAXRBP),
         (RRBP,"Random-RBP",MAXRBP),
         (LRBP,"Local-RBP",MAXRBP),
         (LIST,"List-RBP",MAXRBP)]


############################# FLOODING MODE FLAGS ##############################

# FLOOMODE = "MKAY"
FLOOMODE = "TANH"
# FLOOMODE = "ALTN"
# FLOOMODE = "TABL"
# FLOOMODE = "MSUM"

FAST::Bool = true  # fast flooding update when using tanh mode (default:true)

################################ CONTROL FLAGS #################################

SAVEDATA::Bool = false
PRINTTEST::Bool = false

################################# RBP CONSTANTS ################################

DECAYRBP::Float64 = 0.5
DECAYRRBP::Float64 = 0.7
DECAYLRBP::Float64 = 0.9
DECAYLIST::Float64 = 0.5
decay = Dict(modes[4][2] => DECAYRBP,
             modes[5][2] => DECAYRRBP,
             modes[6][2] => DECAYLRBP,
             modes[7][2] => DECAYLIST)

SAMPLESIZE::Int = 51
LISTSIZE::Int = 200

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
rgn_samples_seeds = zeros(Int,NTHREADS)
for i in eachindex(rgn_noise_seeds)
    rgn_noise_seeds[i] = SEED_NOISE + i - 1
    rgn_samples_seeds[i] = SEED_SAMPL + i - 1
end

#################### JULIA COMPILATION (FOR SPEED) AND TEST ####################
Lr = Dict()
Lq = Dict()
for mode in modes
    if mode[1]
        Lr[mode[2]] , Lq[mode[2]] = performance_simulation(Codeword,SNRTEST,H,
                                        mode[2],1,mode[3],STOP,
                                        rgn_noise_seeds,
                                        printtest=PRINTTEST)
    end
end                             
############################ PERFORMANCE SIMULATION ############################
if TRIALS > 1
    fer_labels = Vector{String}()
    fermax = Vector{Vector{<:AbstractFloat}}()
    FER = Dict()
    BER = Dict()
    for mode in modes
        if mode[1]
            @time FER[mode[2]], BER[mode[2]] = performance_simulation(Codeword,
                                                SNR,H,mode[2],
                                                TRIALS,mode[3],STOP,
                                                rgn_noise_seeds;
                                                rgn_samples_seeds=rgn_samples_seeds)
            push!(fer_labels,mode[2])
            push!(fermax,FER[mode[2]][mode[3],:])
        end            
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

    # FER x Iterations
    labels = Vector{String}()
    for snr in SNR
        push!(labels,"SNR (dB) = $snr")
    end
    labels = permutedims(labels)
    for mode in modes
        if mode[1]
            if mode[2] == "RBP" || mode[2] == "Random-RBP" || mode[2] == "Local-RBP"
                titlefer = "FER $(mode[2]) (decay factor = $(decay[mode[2]]))"
            else
                titlefer = "FER $(mode[2])"
            end
            local p = plot(
                1:mode[3],
                FER[mode[2]],                
                xlabel="Iteration",
                label=labels,
                lw=2,
                title=titlefer,
                ylims=(lim,0)
            )
            display(p)
            SAVEDATA ? savefig(p,"FER_EVOL_"*mode[2]*".png") : nothing
        end
    end
    
    # BER x Iterations
    if !STOP

        for mode in modes
            if mode[1]
                if mode[1] == "RBP" || mode[2] == "Random-RBP" || mode[2] == "Local-RBP"
                    titleber = "BER $(mode[2]) (decay factor = $DECAYRBP)"
                else
                    titleber = "BER $(mode[2])"
                end
                local p = plot(
                    1:mode[3],
                    BER[mode[2]],                
                    xlabel="Iteration",
                    label=labels,
                    lw=2,
                    title=titleber,
                    ylims=(lim-2,0)
                )
                display(p)
                SAVEDATA ? savefig(p,"BER_"*mode[2]*".png") : nothing
            end
        end
    end
end
################################### SAVE DATA ##################################
if TRIALS > 1 && SAVEDATA

    aux = []
    for i in eachindex(FER)
        push!(aux,(fer_labels[i],FER[i]))
    end
    FERS = Dict(aux)
    CSV.write("FERS.csv", DataFrame(FERS), header=true)

    aux = []
    padding = zeros(MAX-MAXRBP)
    for mode in modes
        if mode[1]
            for i in eachindex(SNR)
                title = mode[2] * " (SNR=$(SNR))"
                push!(aux,(title,BER[mode[2]][:,i]))
            end
        end
    end
    BERS = Dict(aux)
    CSV.write("BERS.csv", DataFrame(BERS), header=true)

end