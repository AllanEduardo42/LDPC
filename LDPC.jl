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

############################## 2) INCLUDED FILES ###############################

include("PEG.jl")
include("GF2_functions.jl")
include("IEEE80216e.jl")
include("NR_LDPC_encode.jl")
include("performance_sim.jl")
include("find_girth.jl")

################################# 3) CONSTANTS #################################

const INF = typemax(Int64)
const INFFLOAT = 1e2
const NINFFLOAT = -INFFLOAT
const ALPHA = 0.875               # Min-Sum attenuation factor

now_ = string(now())

# Seeds
SEED_NOISE::Int = 1428
SEED_GRAPH::Int = 5714
SEED_SAMPL::Int = 2857
SEED_MESSA::Int = 1000

############################### 4) CONTROL FLAGS ###############################

TEST::Bool = true
PRIN::Bool = true
MTHR::Bool = true                       
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are 
# not printed)

################################## 5) NUMBERS ##################################

MAX::Int = 50
MAXRBP::Int = 5
DECAY::Float64 = 0.8
SNR = collect(0.8:0.4:1.6)
SNRTEST = [1.6]
TRIALS = 512*[1, 10, 100, 1000, 10000]
TRIALSTEST = [1]

################################ 6) BP SCHEDULE ################################

Modes = zeros(Bool,7)

Modes = ["Flooding","LBP","iLBP","RBP","Local-RBP","List-RBP"]

# BP type: "MKAY", "TANH", "FAST", "ALTN", "TABL", "MSUM"
Bptypes = Vector{String}(undef,7)

# maximum number of BP iterations
Maxiters = zeros(Int,7)

Decays = Vector{Union{Vector{<:AbstractFloat},Nothing}}(nothing,6)

#Flooding
Modes[1] = 0
Bptypes[1] = "FAST"
Maxiters[1] = MAX

#LBP
Modes[2] = 0
Bptypes[2] = "FAST"
Maxiters[2] = MAX

#instantaneos-LBP
Modes[3] = 0
Bptypes[3] = "FAST"
Maxiters[3] = MAX

#RBP
Modes[4] = 0
Bptypes[4] = "FAST"
Maxiters[4] = MAXRBP
Decays[4] = [DECAY]
     
#Local-RBP
Modes[5] = 0
Bptypes[5] = "FAST"
Maxiters[5] = MAXRBP
# Decays[5] = collect(0.0:0.1:1.0)
Decays[5] = [DECAY]

#List-RBP
Modes[6] = 1
Bptypes[6] = "FAST"
Maxiters[6] = MAXRBP
Decays[6] = [DECAY]
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
    for i in eachindex(Modes)
        if Modes[i]
            if Decays[i] !== nothing
                for decay in Decays[i]                
                    name = Modes[i]*" $decay"
                    LR[name] , LQ_[name] = performance_sim(
                        SNRTEST,
                        Modes[i],
                        TRIALSTEST,
                        Maxiters[i],
                        Bptypes[i],
                        decay,
                        STOP,
                        MTHR;
                        test=TEST,
                        printtest = TEST ? PRIN : false)
                end
            else
                LR[Modes[i]] , LQ_[Modes[i]] = performance_sim(
                    SNRTEST,
                    Modes[i],
                    TRIALSTEST,
                    Maxiters[i],
                    Bptypes[i],
                    Decays[i],
                    STOP,
                    MTHR;
                    test=TEST,
                    printtest = TEST ? PRIN : false)
            end
        end
    end                             
else
    Fer_labels = Vector{String}()
    Fermax = Vector{Vector{<:AbstractFloat}}()
    FER = Dict()
    BER = Dict()
    for i in eachindex(Modes)
        if Modes[i]
            if Decays[i] !== nothing
                for decay in Decays[i]                
                    name = Modes[i]*" $decay"
                    FER[name], BER[name] = performance_sim(
                        SNR,
                        Modes[i],
                        TRIALS,
                        Maxiters[i],
                        Bptypes[i],
                        decay,
                        STOP,
                        MTHR)
                    push!(Fer_labels,name*" ($(Bptypes[i]))")
                    push!(Fermax,FER[name][Maxiters[i],:])
                end
            else
                FER[Modes[i]], BER[Modes[i]] = performance_sim(
                    SNR,
                    Modes[i],
                    TRIALS,
                    Maxiters[i],
                    Bptypes[i],
                    Decays[i],
                    STOP,
                    MTHR)
                push!(Fer_labels,Modes[i]*" ($(Bptypes[i]))")
                push!(Fermax,FER[Modes[i]][Maxiters[i],:])
            end
        end            
    end

##################################### SAVE #####################################

    if length(ARGS) == 0
        SAVE = false
    elseif ARGS[1] == "true"
        SAVE = true
    else
        SAVE = false
    end
    
################################### PLOTTING ###################################
    plotlyjs()
    lim = log10(1/maximum(TRIALS))

    # FER x SNR
    if length(SNR) > 1
        Fer_labels = permutedims(Fer_labels)    
        p = plot(
            SNR,Fermax,
            xlabel="SNR (dB)",
            label=Fer_labels,
            lw=2,
            title="FER (Decay factor = $DECAY)",
            ylims=(lim,0)
        )
        SAVE ? savefig(p,"./Saved Data/FER_"*now_*".svg") : display(p) 
    end

    # BER x Iterations
    if !STOP

        # FER x Iterations
        labels = Vector{String}()
        for snr in SNR
            push!(labels,"SNR (dB) = $snr")
        end
        labels = permutedims(labels)
        for i in eachindex(Modes)
            if Modes[i]
                if Decays[i] !== nothing
                    for decay in Decays[i]
                        name = Modes[i]*" $decay"
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
                        SAVE ? savefig(p,"./Saved Data/FER_"*Modes[i]*"_"*now_*".svg") : display(p)
                    end
                else
                    titlefer = "FER $(Modes[i])"                
                    local p = plot(
                        1:Maxiters[i],
                        FER[Modes[i]],                
                        xlabel="Iteration",
                        label=labels,
                        lw=2,
                        title=titlefer,
                        ylims=(lim,0)
                    )
                    SAVE ? savefig(p,"./Saved Data/FER_"*Modes[i]*"_"*now_*".svg") : display(p)
                end
            end
        end

        for i in eachindex(Modes)
            if Modes[i]
                if Decays[i] !== nothing
                    for decay in Decays[i]
                        name = Modes[i]*" $decay"
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
                        SAVE ? savefig(p,"./Saved Data/BER_"*Modes[i]*"_"*now_*".svg") : display(p)
                    end                    
                else
                    titleber = "BER $(Modes[i])"
                    local p = plot(
                            1:Maxiters[i],
                            BER[Modes[i]],                
                            xlabel="Iteration",
                            label=labels,
                            lw=2,
                            title=titleber,
                            ylims=(lim-2,0)
                        )
                        SAVE ? savefig(p,"./Saved Data/BER_"*Modes[i]*"_"*now_*".svg") : display(p)
                end
            end
        end
    end
################################### SAVE DATA ##################################
    if SAVE

        aux = []
        for i in eachindex(Fermax)
            push!(aux,(Fer_labels[i],Fermax[i]))
        end
        FERS = Dict(aux)
        CSV.write("./Saved Data/FERMAX_"*now_*".csv", DataFrame(FERS), header=true)

    end
end