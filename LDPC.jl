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
    now_ = string(now())
    mkdir("./Saved Data/"*now_)
    file = open("./Saved Data/"*now_*"/output.txt", "w")
else
    SAVE = false
end


################################# 3) CONSTANTS #################################

const INF = typemax(Int64)
const INFFLOAT = 1e3
const NINFFLOAT = -INFFLOAT
const CLIP = 0.1*INFFLOAT
const NCLIP = -CLIP
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
DECAY::Float64 = 0.85
SNR = [1.2, 1.6, 1.8, 2.0]
TRIALS = 10 .^(0:length(SNR)-1)*2^10

# TEST
MAXITER_TEST::Int = 1
SNR_TEST = 2.0
TRIALS_TEST = 1
DECAY_TEST = DECAY

################################ 6) BP SCHEDULE ################################

Modes = ["Flooding","LBP","RBP","Local-RBP","List-RBP","Mod-List-RBP",
         "Random-List-RBP"]
Num_modes = length(Modes)
Active = zeros(Bool,Num_modes)
Listsizes = zeros(Int,3)

# BP type: "MKAY", "TANH", "FAST", "ALTN", "TABL", "MSUM"
Bptypes = Vector{String}(undef,Num_modes)

# maximum number of BP iterations
Maxiters = zeros(Int,Num_modes)

Decays = Vector{Vector{<:AbstractFloat}}(undef,Num_modes)
for i in 1:Num_modes
    Decays[i] = [0.0]
end

# Flooding
Active[1] = 0
Bptypes[1] = "FAST"
Maxiters[1] = MAXITER

# LBP
Active[2] = 0
Bptypes[2] = "FAST"
Maxiters[2] = MAXITER

# RBP
Active[3] = 0
Bptypes[3] = "FAST"
Maxiters[3] = MAXIRBP
# Decays[3] = [DECAY]
Decays[3] = [0.7, 0.8, 0.9, 1.0]

# Local-RBP
Active[4] = 0
Bptypes[4] = "FAST"
Maxiters[4] = MAXIRBP
Decays[4] = [DECAY]

# List-RBP
Active[5] = 1
Bptypes[5] = "FAST"
Maxiters[5] = MAXIRBP
Decays[5] = [0.7, 0.8, 0.9, 1.0]

# Mod-List-RBP
Active[6] = 0
Bptypes[6] = "FAST"
Maxiters[6] = MAXIRBP
Decays[6] = [DECAY]    

# Random-List-RBP
Active[7] = 0
Bptypes[7] = "FAST"
Maxiters[7] = MAXIRBP
Decays[7] = [DECAY]

# List-RBP sizes
Listsizes[1] = 16
Listsizes[2] = 1
Listsizes[3] = 16


########################### 7) MESSAGE AND CODEWORD ############################

# Message (Payload) size
A::Int = 1008
# Rate
R::Float64 = 1/2
# LDPC protocol: 1 = NR-LDPC; 2 = PEG; 3 = IEEE80216e;
LDPC::Int = 3
    densities = 2:11

############################# PARITY-CHECK MATRIX #############################

Msg = zeros(Bool,A)

if LDPC == 1
    Zf = 0
    rv = 0
    Cword, H, E_H, nr_ldpc_data = NR_LDPC_encode(Msg,R,rv)
    M::Int, N::Int = size(H)
    girth = find_girth(H,100000)
else
    nr_ldpc_data = NR_LDPC_DATA(0,0,0,0,0,0,0,0,0,0,"0",0,[0],0,[false])
    L = round(Int,A/R)
    if LDPC == 2
        N::Int = L
        M::Int = L - A
        # Vector of the variable node degrees
        D = rand(Xoshiro(SEED_GRAPH),densities,N-M)
        # Generate Parity-Check Matrix by the PEG algorithm
        H, girth = PEG(D,M,N)
    elseif LDPC == 3
        H,Zf,E_H = IEEE80216e(L,R)
        M::Int,N::Int = size(H)
        L = N
        girth = find_girth(H,100000)
    else
        throw(ArgumentError(
                    lazy"CHECK = $CHECK, but must be 1, 2 or 3"
                ))
    end
end

str = 
"""############################### LDPC parameters ################################

Parity Check Matrix: $M x $N"""

println(str)
if SAVE
    println(file,str)
end

display(sparse(H))

str = """

Graph girth = $girth
"""
println(str)
if SAVE
    println(file,str)
end

# Number of Threads
if !TEST
    NTHREADS = Threads.nthreads()
else
    NTHREADS = 1
end

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

Rgn_noise_seeds = zeros(Int,NTHREADS)
Rgn_samples_seeds = zeros(Int,NTHREADS)
Rgn_message_seeds = zeros(Int,NTHREADS)
for i in eachindex(Rgn_noise_seeds)
    Rgn_noise_seeds[i] = SEED_NOISE + i - 1
    Rgn_samples_seeds[i] = SEED_SAMPL + i - 1
    Rgn_message_seeds[i] = SEED_MESSA + 1 - 1
end

############################ PERFORMANCE SIMULATION ############################
if TEST
    if TEST
        LR = Dict()
        LQ_ = Dict()
        Max_residues = Dict()
        for i in eachindex(Active)
            if Active[i]
                LR[Modes[i]], LQ_[Modes[i]], Max_residues[Modes[i]] = performance_sim(
                                            SNR_TEST,
                                            Modes[i],
                                            TRIALS_TEST,
                                            MAXITER_TEST,
                                            Bptypes[i],
                                            DECAY_TEST)
            end
        end
    end
else
    if STOP
        Fer_labels = Vector{String}()
        Fermax = Vector{Vector{<:AbstractFloat}}()
    end
    FER = Dict()
    BER = Dict()
    for i in eachindex(Active)
        if Active[i]
            for decay in Decays[i]
                if decay != 0.0
                    mode = Modes[i]*" $decay"
                else
                    mode = Modes[i]
                end
                FER[mode], BER[mode] = performance_sim(
                                        SNR,
                                        Modes[i],
                                        TRIALS,
                                        Maxiters[i],
                                        Bptypes[i],
                                        decay)
                if STOP
                    push!(Fer_labels,mode*" ($(Bptypes[i]))")
                    push!(Fermax,FER[mode][Maxiters[i],:])
                end
            end
        end
    end  

    ############################### SAVE DATA ##################################
    if SAVE        
        if STOP
            aux = []
            for i in eachindex(Fermax)
                push!(aux,(Fer_labels[i],Fermax[i]))
            end
            FERS = Dict(aux)
            CSV.write("./Saved Data/"*now_*"/FERMAX.csv", DataFrame(FERS), header=true)
        else
            for i in eachindex(Active)
                if Active[i]
                    for decay in Decays[i]
                        if decay != 0.0
                            mode = Modes[i]*" $decay"
                        else
                            mode = Modes[i]
                        end
                        open("./Saved Data/"*now_*"/FER_"*mode*".txt","w") do io
                            writedlm(io,FER[mode])
                        end
                        open("./Saved Data/"*now_*"/BER_"*mode*".txt","w") do io
                            writedlm(io,BER[mode])
                        end
                    end
                end
            end
        end
    else
################################### PLOTTING ###################################
        plotlyjs()
        lim = log10(1/maximum(TRIALS))
    
        # FER x SNR
        if STOP && length(SNR) > 1
            Fer_labels = permutedims(Fer_labels)
            p = plot(
                SNR,Fermax,
                xlabel="SNR (dB)",
                label=Fer_labels,
                lw=2,
                title="FER MAXITER",
                ylims=(lim,0)
            )
            display(p)
        end
    
        if !STOP
    
            # FER x Iterations
            labels = Vector{String}()
            for snr in SNR
                push!(labels,"SNR (dB) = $snr")
            end
            labels = permutedims(labels)
            for i in eachindex(Active)
                if Active[i]
                    for decay in Decays[i]
                        if decay != 0.0
                            mode = Modes[i]*" $decay"
                        else
                            mode = Modes[i]
                        end
                        titlefer = "FER "*mode
                        local p = plot(
                            1:Maxiters[i],
                            FER[mode],
                            xlabel="Iteration",
                            label=labels,
                            lw=2,
                            title=titlefer,
                            ylims=(lim,0)
                        )
                       display(p)
                    end
                end
            end
    
            for i in eachindex(Active)
                if Active[i]
                    for decay in Decays[i]
                        if decay != 0.0
                            mode = Modes[i]*" $decay"
                        else
                            mode = Modes[i]
                        end
                        titleber = "BER $mode"
                        local p = plot(
                            1:Maxiters[i],
                            BER[mode],
                            xlabel="Iteration",
                            label=labels,
                            lw=2,
                            title=titleber,
                            ylims=(lim-2,0)
                        )
                        display(p)
                    end
                end
            end
        end
    end
end
println("The End!")
if SAVE
    println(file, "The End!")
    close(file)
end