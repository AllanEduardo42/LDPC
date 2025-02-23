################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
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
MTHR::Bool = true
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are
# not printed)

################################## 5) NUMBERS ##################################

MAX::Int = 50
MAXRBP::Int = 20
DECAY::Float64 = 0.8
SNR = collect(1.2:0.4:2.0)
SNRTEST = [2.0]
TRIALS = [1, 10, 1000]*2^10
TRIALS = TRIALS[1:length(SNR)]
TRIALSTEST = [2]

################################ 6) BP SCHEDULE ################################

Modes = zeros(Bool,7)

Names = ["Flooding","LBP","RBP","Local-RBP","List-RBP"]

# BP type: "MKAY", "TANH", "FAST", "ALTN", "TABL", "MSUM"
Bptypes = Vector{String}(undef,6)

# maximum number of BP iterations
Maxiters = zeros(Int,5)

Decays = Vector{Vector{<:AbstractFloat}}(undef,5)
for i in eachindex(Decays)
    Decays[i] = [0.0]
end

#Flooding
Modes[1] = 0
Bptypes[1] = "FAST"
Maxiters[1] = MAX

#LBP
Modes[2] = 0
Bptypes[2] = "FAST"
Maxiters[2] = MAX

#RBP
Modes[3] = 1
Bptypes[3] = "FAST"
Maxiters[3] = MAXRBP
Decays[3] = [DECAY]
# Decays[4] = collect(1.0:-0.1:0.8)

#Local-RBP
Modes[4] = 0
Bptypes[4] = "FAST"
Maxiters[4] = MAXRBP
Decays[4] = [DECAY]

#List-RBP
Modes[5] = 0
Bptypes[5] = "FAST"
Maxiters[5] = MAXRBP
Decays[5] = [DECAY]
    # List-RBP size
    LISTSIZE::UInt = 16
    LISTSIZE2::UInt = 1

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
if MTHR && !TEST
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
    LR = Dict()
    LQ_ = Dict()
    Max_residues = Dict()
    for i in eachindex(Modes)
        if Modes[i]
            for decay in Decays[i]
                if decay != 0.0
                    name = Names[i]*" $decay"
                else
                    name = Names[i]
                end
                LR[name], LQ_[name], Max_residues[name] = performance_sim(
                                            SNRTEST,
                                            Names[i],
                                            TRIALSTEST,
                                            Maxiters[i],
                                            Bptypes[i],
                                            decay,
                                            STOP,
                                            MTHR;
                                            test=TEST,
                                            printtest = TEST ? PRIN : false)
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
    for i in eachindex(Modes)
        if Modes[i]
            for decay in Decays[i]
                if decay != 0.0
                    name = Names[i]*" $decay"
                else
                    name = Names[i]
                end
                FER[name], BER[name] = performance_sim(
                    SNR,
                    Names[i],
                    TRIALS,
                    Maxiters[i],
                    Bptypes[i],
                    decay,
                    STOP,
                    MTHR)
                if STOP
                    push!(Fer_labels,name*" ($(Bptypes[i]))")
                    push!(Fermax,FER[name][Maxiters[i],:])
                end
            end
        end
    end

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
            title="FER MAX",
            ylims=(lim,0)
        )
        SAVE ? savefig(p,"./Saved Data/"*now_*"/FERMAX.svg") : display(p)
    end

    if !STOP

        # FER x Iterations
        labels = Vector{String}()
        for snr in SNR
            push!(labels,"SNR (dB) = $snr")
        end
        labels = permutedims(labels)
        for i in eachindex(Modes)
            if Modes[i]
                for decay in Decays[i]
                    if decay != 0.0
                        name = Names[i]*" $decay"
                    else
                        name = Names[i]
                    end
                    titlefer = "FER "*name
                    local p = plot(
                        1:Maxiters[i],
                        FER[name],
                        xlabel="Iteration",
                        label=labels,
                        lw=2,
                        title=titlefer,
                        ylims=(lim,0)
                    )
                    SAVE ? savefig(p,"./Saved Data/"*now_*"/FER_"*name*".svg") : display(p)
                end
            end
        end

        for i in eachindex(Modes)
            if Modes[i]
                for decay in Decays[i]
                    if decay != 0.0
                        name = Names[i]*" $decay"
                    else
                        name = Names[i]
                    end
                    titleber = "BER $name"
                    local p = plot(
                        1:Maxiters[i],
                        BER[name],
                        xlabel="Iteration",
                        label=labels,
                        lw=2,
                        title=titleber,
                        ylims=(lim-2,0)
                    )
                    SAVE ? savefig(p,"./Saved Data/"*now_*"/BER_"*name*".svg") : display(p)
                end
            end
        end
    end
################################### SAVE DATA ##################################
    if SAVE        
        if STOP
            aux = []
            for i in eachindex(Fermax)
                push!(aux,(Fer_labels[i],Fermax[i]))
            end
            FERS = Dict(aux)
            CSV.write("./Saved Data/"*now_*"/FERMAX.csv", DataFrame(FERS), header=true)
        else
            for i in eachindex(Modes)
                if Modes[i]
                    for decay in Decays[i]
                        if decay != 0.0
                            name = Names[i]*" $decay"
                        else
                            name = Names[i]
                        end
                        open("./Saved Data/"*now_*"/FER_"*name*".txt","w") do io
                            writedlm(io,FER[name])
                        end
                        open("./Saved Data/"*now_*"/BER_"*name*".txt","w") do io
                            writedlm(io,BER[name])
                        end
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