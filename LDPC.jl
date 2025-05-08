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
const TABLESIZE = 8192
const TABLERANGE = 10
const SIZE_PER_RANGE = TABLESIZE/TABLERANGE

# Seeds
SEED_NOISE::Int = 1428
SEED_GRAPH::Int = 5714
SEED_MESSA::Int = 1000

############################### 4) CONTROL FLAGS ###############################

TEST::Bool = true
PRIN::Bool = true
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are
# not printed)

################################## 5) NUMBERS ##################################

MAXITER::Int = 20
# FACTORS = [0.7, 0.8, 0.9, 1.0]
# FACTORS = collect(0.1:0.1:1.0)
FACTORS = [1.0]
# EbN0 = [1.2, 1.4, 1.6, 1.8]
EbN0 = [1.5]
TRIALS = 10 .^(0:length(EbN0)-1)*2^14
RELATIVE::Bool = false

# TEST
MAXITER_TEST::Int = 1
EbN0_TEST::Float64 = 1.5
TRIALS_TEST::Int = 1
DECAY_TEST::Float64 = 1.0

################################ 6) BP SCHEDULE ################################

MODES = ["Flooding","LBP","RBP","List-RBP","VN-RBP","Genius-RBP","NW-RBP"]
NUM_MODES = length(MODES)
ACTIVE = zeros(Bool,NUM_MODES)
LISTSIZES = zeros(Int,4)

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
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# LBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "TANH"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Variable Node RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Genius-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# NW-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP sizes (min values = 4 and 2)
LISTSIZES[1] = 128
LISTSIZES[2] = 16

############################### 7) LOOKUP TABLE ################################

PHI = lookupTable()

########################### 8) MESSAGE AND CODEWORD ############################
# Rate
RR::Float64 = 1/2
# Message (Payload) size
# AA::Int = 576*RR
AA::Int = 1008
# LDPC protocol: NR5G = NR-LDPC (5G); PEG = PEG; WiMAX = IEEE80216e;
PROTOCOL::String = "NR5G"
    DENSITIES = 1:8

############################# PARITY-CHECK MATRIX #############################

MSG = zeros(Bool,AA)

if PROTOCOL == "NR5G"
    ZF = 0
    RV = 0
    CWORD, HH, E_H, NR_LDPC_DATA = NR_LDPC_encode(MSG,RR,RV)
    MM, NN = size(HH)
    GIRTH = find_girth(HH,100000)
    GG = nothing
else
    NR_LDPC_DATA = nothing
    LL = round(Int,AA/RR)
    if PROTOCOL == "PEG"
        ZF = 0
        E_H = nothing
        NN = LL
        MM = LL - AA
        # Vector of the variable node degrees
        DD = rand(Xoshiro(SEED_GRAPH),DENSITIES,NN)
        # Generate Parity-Check Matrix by the PEG algorithm
        if AA == 1008 && RR == 1/2
            HH = readdlm("./PEG_1008_2016.txt",'\t', Bool,'\n')
            GIRTH = find_girth(HH,100000)
        else
            HH, GIRTH = PEG(DD,MM,NN)
        end
        KK = NN-MM
        GG = [I(KK); gf2_mat_mult(gf2_inverse(HH[:,KK+1:NN]), HH[:,1:KK])]
    elseif PROTOCOL == "WiMAX"
        # N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
        # 1632,1728,1824,1920,2016,2112,2208,2304}.    
        # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.
        HH,ZF,E_H = IEEE80216e(LL,RR,"A")
        MM,NN = size(HH)
        LL = NN
        GIRTH = find_girth(HH,100000)
        GG = nothing
    end
end

# list of checks and variables nodes
NC  = make_cn2vn_list(HH)
NV  = make_vn2cn_list(HH)

STR = 
"""############################### LDPC parameters ################################
LDPC Protocol: """
if PROTOCOL == "NR5G"
    STR *= "NR-LDPC (5G)"
elseif PROTOCOL == "PEG"
    STR *= "PEG"
elseif PROTOCOL == "WiMAX"
    STR *= "IEEE80216e (WiMAX)"
end
STR *= "\nParity Check Matrix: $MM x $NN"
println(STR)
if SAVE
    println(FILE,STR)
end

display(sparse(HH))

STR = """

Graph girth = $GIRTH
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
RGN_MESSAGE_SEEDS = zeros(Int,NTHREADS)
for i in eachindex(RGN_NOISE_SEEDS)
    RGN_NOISE_SEEDS[i] = SEED_NOISE + i - 1
    RGN_MESSAGE_SEEDS[i] = SEED_MESSA + 1 - 1
end

############################ PERFORMANCE SIMULATION ############################
if TEST
    if TEST
        LRM = Dict()
        LQM = Dict()
        for i in eachindex(ACTIVE)
            if ACTIVE[i]
                LRM[MODES[i]], LQM[MODES[i]] = performance_sim(
                                            EbN0_TEST,
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
                                        EbN0,
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
        # gr()
        LIMFER = 1/maximum(TRIALS)
        LIMBER = 1/(maximum(TRIALS)*NN)
    
        # FER x EbN0
        if STOP && length(EbN0) > 1
            FER_LABELS = permutedims(FER_LABELS)
            PLOT = plot(
                EbN0,FERMAX,
                xlabel="EbN0 (dB)",
                label=FER_LABELS,
                lw=2,
                title="FER MAXITER",
                ylims=(LIMFER,0)
            )
            display(PLOT)
        end
    
        if !STOP  
            FB = ["F","B"]  
            # FER x Iterations
            for j=1:2
                p = plot()
                title = FB[j]*"ER $PROTOCOL (rate = $RR)"
                for i in eachindex(ACTIVE)
                    if ACTIVE[i]
                        for decay in DECAYS[i]
                            if decay != 0.0
                                mode = MODES[i]*" $decay"
                            else
                                mode = MODES[i]
                            end
                            labels = Vector{String}()
                            for enr in EbN0
                                push!(labels,"$mode, Eb/N0 = $(enr)dB")
                            end
                            labels = permutedims(labels)
                            if FB[j] == "F"
                                y = log10.(FER[mode])
                                lim = log10(LIMFER)
                            else
                                y = log10.(BER[mode])
                                lim = log10(LIMBER)
                            end
                            x = 1:MAXITERS[i]
                            p = plot!(
                                x,
                                y,
                                # yscale=:log10,
                                xlabel="Iteration",
                                # minorgrid=true,
                                label=labels,
                                lw=2,
                                title=title,
                                ylims=(lim,0)
                            )
                        end
                    end
                end
                display(p)
            end
            # PLOT = plot!(1:MAXITERS[i],
            #       ones(Int,MAXITER)*LIMFER,
            #       label = "minimum",
            #       lw=2,
            #       color=:black,
            #       ls=:dash)
        end
    end
end
println("The End!")
if SAVE
    println(FILE, "The End!")
    close(FILE)
end