################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# LDPC coding-decoding performance simulation using Believe Propagation (BP)

############################## 1) JULIA PACKAGES ###############################

using LinearAlgebra
using Random
using Statistics
using Plots
using SparseArrays
using CSV, DataFrames
using Dates
using DelimitedFiles

############################## 2) INCLUDED FILES ###############################

include("PEG.jl")
include("GF2_functions.jl")
include("IEEE80216e.jl")
include("NR_LDPC_encode.jl")
include("performance_sim.jl")
include("find_girth.jl")

##################################### SAVE #####################################

if length(ARGS) == 0
    SAVE = false
elseif ARGS[1] == "true"
    SAVE = true
    NOW = string(now())
    mkdir("./Saved Data/"*NOW)
    FILE = open("./Saved Data/"*NOW*"/output.txt", "w")
else
    SAVE = false
end


################################# 3) CONSTANTS #################################

const INF = typemax(Int64)
const INFFLOAT = 1e3
const NINFFLOAT = -INFFLOAT
const ALPHA = 0.875               # Min-Sum attenuation factor

# Seeds
SEED_NOISE::Int = 1428
SEED_GRAPH::Int = 5714
SEED_SAMPL::Int = 2857
SEED_MESSA::Int = 1000

############################### 4) CONTROL FLAGS ###############################

TEST::Bool = false
PRIN::Bool = true
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are
# not printed)

################################## 5) NUMBERS ##################################

MAXITER::Int = 50
MAXIRBP::Int = 30
# FACTORS = [0.7, 0.8, 0.9, 1.0]
FACTORS = [0.9]
SNR = [1.2, 1.6]
TRIALS = 10 .^(0:length(SNR)-1)*2^7

# TEST
MAXITER_TEST::Int = 1
SNR_TEST::Float64 = 2.0
TRIALS_TEST::Int = 1
DECAY_TEST::Float64 = 0.9

################################ 6) BP SCHEDULE ################################

MODES = ["Flooding","LBP","VLBP","RBP","Local-RBP","List-RBP","Mod-List-RBP",
         "Random-List-RBP","List-RBP-genius","NRBP"]
NUM_MODES = length(MODES)
ACTIVE = zeros(Bool,NUM_MODES)
LISTSIZES = zeros(Int,3)

# BP type: "MKAY", "TANH", "FAST", "ALTN", "TABL", "MSUM"
BPTYPES = Vector{String}(undef,NUM_MODES)

# maximum number of BP iterations
MAXITERS = zeros(Int,NUM_MODES)

DECAYS = Vector{Vector{<:AbstractFloat}}(undef,NUM_MODES)
for i in 1:NUM_MODES
    DECAYS[i] = [0.0]
end

i = 1
# Flooding
ACTIVE[i] = 1
BPTYPES[i] = "MKAY"
MAXITERS[i] = MAXITER

# LBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# VLBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Local-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Mod-List-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Random-List-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP-genius
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Node-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = [1.0]

# List-RBP sizes
LISTSIZES[1] = 16
LISTSIZES[2] = 2
LISTSIZES[3] = 16


########################### 7) MESSAGE AND CODEWORD ############################

# Message (Payload) size
AA::Int = 1008
# Rate
RR::Float64 = 1/2
# LDPC protocol: 1 = NR-LDPC; 2 = PEG; 3 = IEEE80216e;
LDPC::Int = 1
    DENSITIES = 2:12

############################# PARITY-CHECK MATRIX #############################

MSG = zeros(Bool,AA)

if LDPC == 1
    ZF = 0
    RV = 0
    CWORD, HH, E_H, NR_LDPC_DATA = NR_LDPC_encode(MSG,RR,RV)
    MM, NN = size(HH)
    GIRTH = find_girth(HH,100000)
    GG = nothing
else
    NR_LDPC_DATA = nothing
    LL = round(Int,AA/RR)
    if LDPC == 2
        E_H = nothing
        NN = LL
        MM = LL - AA
        # Vector of the variable node degrees
        DD = rand(Xoshiro(SEED_GRAPH),DENSITIES,NN)
        # Generate Parity-Check Matrix by the PEG algorithm
        @time HH, GIRTH = PEG(DD,MM,NN)
        GG = [gf2_mat_mult(gf2_inverse(HH[:,1:MM]),HH[:,MM+1:NN]); I(NN-MM)]
    elseif LDPC == 3
        # N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
        # 1632,1728,1824,1920,2016,2112,2208,2304}.    
        # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.
        HH,ZF,E_H = IEEE80216e(LL,RR)
        MM,NN = size(HH)
        LL = NN
        GIRTH = find_girth(HH,100000)
        GG = nothing
    end
end

# list of checks and variables nodes
CN2VN  = make_cn2vn_list(HH)
VN2CN  = make_vn2cn_list(HH)

STR = 
"""############################### LDPC parameters ################################

Parity Check Matrix: $MM x $NN"""

println(STR)
if SAVE
    println(FILE,STR)
end

display(sparse(HH))

STR = """

Graph GIRTH = $GIRTH
"""
println(STR)
if SAVE
    println(FILE,STR)
end

# Number of Threads
if !TEST
    NTHREADS = Threads.nthreads()
else
    NTHREADS = 1
end

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

RGN_NOISE_SEEDS = zeros(Int,NTHREADS)
RGN_SAMPLE_SEEDS = zeros(Int,NTHREADS)
RGN_MESSAGE_SEEDS = zeros(Int,NTHREADS)
for i in eachindex(RGN_NOISE_SEEDS)
    RGN_NOISE_SEEDS[i] = SEED_NOISE + i - 1
    RGN_SAMPLE_SEEDS[i] = SEED_SAMPL + i - 1
    RGN_MESSAGE_SEEDS[i] = SEED_MESSA + 1 - 1
end

############################ PERFORMANCE SIMULATION ############################
if TEST
    if TEST
        LRM = Dict()
        LQM = Dict()
        MAX_RESIDUES = Dict()
        for i in eachindex(ACTIVE)
            if ACTIVE[i]
                LRM[MODES[i]], LQM[MODES[i]], MAX_RESIDUES[MODES[i]] = performance_sim(
                                            SNR_TEST,
                                            MODES[i],
                                            TRIALS_TEST,
                                            MAXITER_TEST,
                                            BPTYPES[i],
                                            DECAY_TEST)
            end
        end
    end
else
    if STOP
        FER_LABELS = Vector{String}()
        FERMAX = Vector{Vector{<:AbstractFloat}}()
    end
    FER = Dict()
    BER = Dict()
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            for decay in DECAYS[i]
                if decay != 0.0
                    mode = MODES[i]*" $decay"
                else
                    mode = MODES[i]
                end
                FER[mode], BER[mode] = performance_sim(
                                        SNR,
                                        MODES[i],
                                        TRIALS,
                                        MAXITERS[i],
                                        BPTYPES[i],
                                        decay)
                if STOP
                    push!(FER_LABELS,mode*" ($(BPTYPES[i]))")
                    push!(FERMAX,FER[mode][MAXITERS[i],:])
                end
            end
        end
    end  

    ############################### SAVE DATA ##################################
    if SAVE        
        if STOP
            aux = []
            for i in eachindex(FERMAX)
                push!(aux,(FER_LABELS[i],FERMAX[i]))
            end
            FERS = Dict(aux)
            CSV.write("./Saved Data/"*NOW*"/FERMAX.csv", DataFrame(FERS), header=true)
        else
            for i in eachindex(ACTIVE)
                if ACTIVE[i]
                    for decay in DECAYS[i]
                        if decay != 0.0
                            mode = MODES[i]*" $decay"
                        else
                            mode = MODES[i]
                        end
                        open("./Saved Data/"*NOW*"/FER_"*mode*".txt","w") do io
                            writedlm(io,FER[mode])
                        end
                        open("./Saved Data/"*NOW*"/BER_"*mode*".txt","w") do io
                            writedlm(io,BER[mode])
                        end
                    end
                end
            end
        end
    else
################################### PLOTTING ###################################
        plotlyjs()
        LIM = log10(1/maximum(TRIALS))
    
        # FER x SNR
        if STOP && length(SNR) > 1
            FER_LABELS = permutedims(FER_LABELS)
            PLOT = plot(
                SNR,FERMAX,
                xlabel="SNR (dB)",
                label=FER_LABELS,
                lw=2,
                title="FER MAXITER",
                ylims=(LIM,0)
            )
            display(PLOT)
        end
    
        if !STOP    
            # FER x Iterations
            PLOT = plot()
            titlefer = "FER"
            for i in eachindex(ACTIVE)
                if ACTIVE[i]
                    for decay in DECAYS[i]
                        if decay != 0.0
                            mode = MODES[i]*" $decay"
                        else
                            mode = MODES[i]
                        end
                        LABELS = Vector{String}()
                        for snr in SNR
                            push!(LABELS,"$mode, SNR (dB) = $snr")
                        end
                        LABELS = permutedims(LABELS)
                        local PLOT = plot!(
                            1:MAXITERS[i],
                            FER[mode],
                            xlabel="Iteration",
                            label=LABELS,
                            lw=2,
                            title=titlefer,
                            ylims=(LIM,0)
                        )
                    end
                end
            end
            display(PLOT)
            PLOT = plot()
            titleber = "BER"
            for i in eachindex(ACTIVE)
                if ACTIVE[i]
                    for decay in DECAYS[i]
                        if decay != 0.0
                            mode = MODES[i]*" $decay"
                        else
                            mode = MODES[i]
                        end
                        LABELS = Vector{String}()
                        for snr in SNR
                            push!(LABELS,"$mode, SNR (dB) = $snr")
                        end
                        LABELS = permutedims(LABELS)
                        local PLOT = plot!(
                            1:MAXITERS[i],
                            BER[mode],
                            xlabel="Iteration",
                            label=LABELS,
                            lw=2,
                            title=titleber,
                            ylims=(LIM-2,0)
                        )
                    end
                end
            end
            display(PLOT)
        end
    end
end
println("The End!")
if SAVE
    println(FILE, "The End!")
    close(FILE)
end