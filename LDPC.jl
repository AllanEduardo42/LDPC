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
include("NR LDPC/NR_LDPC_parameters.jl")
include("performance_sim.jl")
include("find_girth.jl")
include("LU_encoding.jl")

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
const ALPHA = 0.7              # Min-Sum attenuation factor
const TABLESIZE = 8192
const TABLERANGE = 10
const SIZE_PER_RANGE = TABLESIZE/TABLERANGE

# Seeds
SEED_NOISE::Int = 1428
SEED_MESSA::Int = 1000

############################### 4) CONTROL FLAGS ###############################

TEST::Bool = false
PRIN::Bool = true
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are
# not printed)

################################## 5) NUMBERS ##################################

MAXITER::Int = 10
# FACTORS = [0.7, 0.8, 0.9, 1.0]
# FACTORS = collect(0.1:0.1:1.0)
FACTORS = [0.9]
# EbN0 = [1.2, 1.4, 1.6, 1.8]
EbN0 = [2.5]
TRIALS = 10 .^(0:length(EbN0)-1)*2^14
RELATIVE::Bool = true

# TEST
MAXITER_TEST::Int = 1
EbN0_TEST::Float64 = 1.0
TRIALS_TEST::Int = 10
DECAY_TEST::Float64 = 1.0

################################ 6) BP SCHEDULE ################################

MODES = ["Flooding","LBP","RBP","List-RBP","NW-RBP","VN-RBP"]
NUM_MODES = length(MODES)
ACTIVE = zeros(Bool,NUM_MODES)
LISTSIZES = zeros(Int,4)

# BP type: "MKAY", "TANH", "FAST", "TABL", "MSUM"
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
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# LBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# NW-RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Variable Node RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP sizes (min values = 4 and 2)
LISTSIZES[1] = 16
LISTSIZES[2] = 2

########################### 8) MESSAGE AND CODEWORD ############################

# Message (Payload) size
GG = 576
# Effective Rate
RR = 1/2 - 16/GG  # WiMAX compatibility offset
# LDPC protocol: NR5G = NR-LDPC (5G); PEG = PEG; WiMAX = IEEE80216e;
PROTOCOL::String = "PEG"
    LAMBDA = [0.21, 0.25, 0.25, 0.29, 0]
    RO = [1.0, 0, 0, 0, 0, 0]

############################# PARITY-CHECK MATRIX #############################

if PROTOCOL == "NR5G"
    RV = 0
    AA, KK, RR, G_CRC, LIFTSIZE, NR_LDPC_DATA = NR_LDPC_parameters(GG,RR,RV,false)
    HH, E_H = NR_LDPC_make_parity_check_matrix(LIFTSIZE,
                                               NR_LDPC_DATA.iLS,
                                               NR_LDPC_DATA.bg,
                                               NR_LDPC_DATA.P,
                                               NR_LDPC_DATA.K_prime,
                                               NR_LDPC_DATA.K,
                                               NR_LDPC_DATA.P_Zc)
    MM, NN = size(HH)
    GIRTH = find_girth(HH,100000)
    LL = nothing
    UU = nothing
else
    NN = GG
    AA = round(Int,GG*RR)
    KK, g_CRC = get_CRC_poly(AA)
    _, G_CRC = get_CRC_poly(AA)  
    if PROTOCOL == "PEG"     
        LIFTSIZE = 0
        E_H = nothing
        MM = NN - KK
        # Generate Parity-Check Matrix by the PEG algorithm
        H_PEG, GIRTH = PEG(LAMBDA,RO,MM,NN)
        HH, LL, UU = remake_H(H_PEG,0)
    elseif PROTOCOL == "WiMAX"
        # N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
        # 1632,1728,1824,1920,2016,2112,2208,2304}.    
        # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.
        HH, LIFTSIZE, E_H = IEEE80216e(NN,RR,"A")
        MM,_ = size(HH)
        KK = NN - MM
        AA = KK - 16 # since max(AA) = 2304*5/6 ≤ 3824, g_cRC = {CRC16}      
        GIRTH = find_girth(HH,100000)
        LL = nothing
        UU = nothing
    end
end

# list of checks and variables nodes
NC = make_cn2vn_list(HH)
NV = make_vn2cn_list(HH)

STR = 
"""############################### LDPC parameters ################################
LDPC Protocol: """
if PROTOCOL == "NR5G"
    STR *= "NR-LDPC (5G), Zc = $(NR_LDPC_DATA.Zc)"
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
Effective rate = $(round(RR,digits=3))
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
    LRM = zeros(MM,NN)
    LQM = zeros(MM,NN)
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            global LRM, LQM = performance_sim(EbN0_TEST,
                                              MODES[i],
                                              TRIALS_TEST,
                                              MAXITER_TEST,
                                              BPTYPES[i],
                                              DECAY_TEST)
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
                title = FB[j]*"ER $PROTOCOL (R = $(round(RR,digits=2)))"
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