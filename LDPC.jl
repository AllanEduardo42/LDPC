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
const ALPHA2 = 2*ALPHA 


name = string(now())

# Seeds
SEED_NOISE::Int = 1428
SEED_GRAPH::Int = 5714
SEED_SAMPL::Int = 2857
SEED_MESSA::Int = 1000

############################### 4) CONTROL FLAGS ###############################

MTHR::Bool = true                       
SAVE::Bool = false
PRIN::Bool = true
STOP::Bool = true # stop simulation at zero syndrome (if true, BER curves are 
# not printed)

################################## 5) NUMBERS ##################################

TRIALS::Int = 102400
MAX::Int = 50
MAXRBP::Int = 50
DECAY::Float64 = 0.9
SNRTEST = [3]
SNR = collect(1:0.4:2.2)

################################ 6) BP SCHEDULE ################################

Modes = zeros(Bool,7)

Names = ["Flooding","LBP","iLBP","RBP","Local-RBP","List-RBP"]

# BP type: "MKAY", "TANH", "FAST", "ALTN", "TABL", "MSUM"
Bptypes = Vector{String}(undef,7)

# maximum number of BP iterations
Maxiters = zeros(Int,7)

Decays = zeros(7)

#Flooding
Modes[1] = 1
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
Modes[4] = 1
Bptypes[4] = "FAST"
Maxiters[4] = MAXRBP
Decays[4] = DECAY
     
#Local-RBP
Modes[5] = 0
Bptypes[5] = "MSUM"
Maxiters[5] = MAXRBP
Decays[5] = DECAY

#List-RBP
Modes[6] = 0
Bptypes[6] = "FAST"
Maxiters[6] = MAXRBP
Decays[6] = DECAY
    # List-RBP size
    LISTSIZE::UInt = 64
    LISTSIZE2::UInt = 4

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
    rv = 0
    H, Cword, Zc = NR_LDPC_encode(Msg,R,rv)
    M::Int, N::Int = size(H)
    girth = find_girth(H,100000)
else
    Zc = 0
    L = round(Int,A/R)
    if LDPC == 2
        N::Int = L           
        M::Int = L - A
        # Vector of the variable node degrees
        D = rand(Xoshiro(SEED_GRAPH),densities,N-M)
        # Generate Parity-Check Matrix by the PEG algorithm
        H, girth = PEG(D,M,N)
    elseif LDPC == 3
        H,zf,E_H = IEEE80216e(L,R)
        M::Int,N::Int = size(H)
        L = N
        girth = find_girth(H,100000)
    else
        throw(ArgumentError(
                    lazy"CHECK = $CHECK, but must be 1, 2 or 3"
                ))
    end
end

# Number of Threads
NTHREADS::Int = min(Threads.nthreads(),TRIALS)

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

Rgn_noise_seeds = zeros(Int,NTHREADS)
Rgn_samples_seeds = zeros(Int,NTHREADS)
Rgn_message_seeds = zeros(Int,NTHREADS)
for i in eachindex(Rgn_noise_seeds)
    Rgn_noise_seeds[i] = SEED_NOISE + i - 1
    Rgn_samples_seeds[i] = SEED_SAMPL + i - 1
    Rgn_message_seeds[i] = SEED_MESSA + 1 - 1
end

#################### JULIA COMPILATION (FOR SPEED) AND TEST ####################
LR = Dict()
LQ = Dict()
p = (TRIALS ≤ 2) ? PRIN : false
for i in eachindex(Modes)
    if Modes[i]
        LR[Names[i]] , LQ[Names[i]] = performance_sim(
            SNRTEST,
            Names[i],
            min(TRIALS,2),
            Maxiters[i],
            Bptypes[i],
            Decays[i];
            printtest=p)
    end
end                             
############################ PERFORMANCE SIMULATION ############################
if TRIALS > 2
    Fer_labels = Vector{String}()
    Fermax = Vector{Vector{<:AbstractFloat}}()
    FER = Dict()
    BER = Dict()
    for i in eachindex(Modes)
        if Modes[i]
            @time FER[Names[i]], BER[Names[i]] = performance_sim(
                SNR,
                Names[i],
                TRIALS,
                Maxiters[i],
                Bptypes[i],
                Decays[i])

            push!(Fer_labels,Names[i]*" ($(Bptypes[i]))")
            push!(Fermax,FER[Names[i]][Maxiters[i],:])
        end            
    end
end
################################### PLOTTING ###################################
if TRIALS > 2
    plotlyjs()
    lim = log10(1/TRIALS)

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
        SAVE ? savefig(p,"./Saved Data/FER_"*name*".svg") : display(p) 
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
                if Names[i] == "RBP" || Names[i] == "Random-RBP" ||
                Names[i] == "Local-RBP" || Names[i] == "List-RBP"
                    titlefer = "FER $(Names[i]) (decay factor = $(decay[Names[i]]))"
                else
                    titlefer = "FER $(Names[i])"
                end
                local p = plot(
                    1:Maxiters[i],
                    FER[Names[i]],                
                    xlabel="Iteration",
                    label=labels,
                    lw=2,
                    title=titlefer,
                    ylims=(lim,0)
                )
                SAVE ? savefig(p,"./Saved Data/FER_"*Names[i]*"_"*name*".svg") : display(p) 
            end
        end

        for i in eachindex(Modes)
            if Modes[i]
                if Names[i] == "RBP" || Names[i] == "Random-RBP" || 
                   Names[i] == "Local-RBP" || Names[i] == "List-RBP"
                    titleber = "BER $(Names[i]) (decay factor = $(decay[Names[i]]))"
                else
                    titleber = "BER $(Names[i])"
                end
                local p = plot(
                    1:Maxiters[i],
                    BER[Names[i]],                
                    xlabel="Iteration",
                    label=labels,
                    lw=2,
                    title=titleber,
                    ylims=(lim-2,0)
                )
                SAVE ? savefig(p,"./Saved Data/BER_"*Names[i]*"_"*name*".svg") : display(p)
            end
        end
    end
end
################################### SAVE DATA ##################################
if TRIALS > 2 && SAVE

    aux = []
    for i in eachindex(Fermax)
        push!(aux,(Fer_labels[i],Fermax[i]))
    end
    FERS = Dict(aux)
    CSV.write("./Saved Data/FERMAX_"*name*".csv", DataFrame(FERS), header=true)

end