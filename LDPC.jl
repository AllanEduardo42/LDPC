################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# LDPC coding-decoding performance simulation

################################ JULIA PACKAGES ################################

using LinearAlgebra
using Random
using Statistics
using Plots
using SparseArrays

################################ BP MODE FLAGS ################################

 LBP = false
iLBP = true
 RBP = false
LRBP = false

############################# FLOODING MODE FLAGS ##############################
MKAY = false
TANH = false
ALTN = false
TABL = false
MSUM = false

# fast flooding update when using tanh mode (default:true)
FAST = true

########################## PRINTING IN TEST MODE FLAG ##########################
PRINTING = true

####################### STOP WHEN SYNDROME IS ZERO FLAG ########################
STOP = true

################################### PLOTTING ###################################

PLOT_BER = true
HISTOGRAMS = false

################################ INCLUDED FILES ################################

include("performance_simulation.jl")
include("PEG.jl")
include("GF2_functions.jl")

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1428
SEED2::Int64 = 5714

SIZE::Int64 = 1024
RANGE::Int64 = 20

SIZE_per_RANGE::Float64 = SIZE/RANGE

NREALS::Int = 1
MAX::Int = 30
MAX_RBP::Int = 5

LR_idx::Int = 4;
DECAYCTE::Float64 = 0.9
SAMPLESIZE::Int = 52

#################################### NOISE #####################################

SNR = collect(0:1:8)

############################# PARITY-CHECK MATRIX #############################

### PEG Algorithm

N::Int64 = 512
M::Int64 = 256
Random.seed!(SEED2)

D = rand([2,3,4],N)

H, girth = PEG!(D,M)

G = gf2_nullspace(H)

############################# MESSAGE AND CODEWORD #############################

Message = rand(Bool,N-M)

Codeword = gf2_mat_mult(Matrix(G), Message)

######################### PRINT INFORMATION ON SCREEN #########################

println()
println(
"############################### LDPC parameters ###############################"
)
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

############################## JULIA COMPILATION ###############################
if MKAY
    Lr_mkay, Lq_mkay = performance_simulation(Codeword,[SNR[LR_idx]],H,"MKAY",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if TANH
    Lr_tanh, Lq_tanh = performance_simulation(Codeword,[SNR[LR_idx]],H,"TANH",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if ALTN
    Lr_altn, Lq_altn = performance_simulation(Codeword,[SNR[LR_idx]],H,"ALTN",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if TABL
    Lr_tabl, Lq_tabl = performance_simulation(Codeword,[SNR[LR_idx]],H,"TABL",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if MSUM
    Lr_msum, Lq_msum = performance_simulation(Codeword,[SNR[LR_idx]],H,"MSUM",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if  LBP
      Lr_lbp, Lq_lbp = performance_simulation(Codeword,[SNR[LR_idx]],H, "LBP",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if iLBP
    Lr_ilbp, Lq_ilbp = performance_simulation(Codeword,[SNR[LR_idx]],H,"iLBP",1,MAX;
                                              printing=PRINTING, stop=STOP)
end
if RBP
    Lr_rbp, Lq_rbp, Edges = 
                       performance_simulation(Codeword,[SNR[LR_idx]],H, "RBP",1,MAX_RBP;
                                              printing=PRINTING, stop=STOP)
end
if LRBP
    Lr_lrbp, Lq_lrbp, lEdges = 
                       performance_simulation(Codeword,[SNR[LR_idx]],H,"LRBP",1,MAX_RBP;
                                              printing=PRINTING, stop=STOP)
end
                             
############################ PERFORMANCE SIMULATION ############################
if NREALS > 1
    if MKAY
        @time FER, BER, Iters = 
            performance_simulation(Codeword,SNR,H,"MKAY",NREALS,MAX)
    end
    if TANH
        @time FER_tanh, BER_tanh, Iters_tanh = 
            performance_simulation(Codeword,SNR,H,"TANH",NREALS,MAX)
    end
    if ALTN
        @time FER_altn, BER_altn, Iters_altn = 
            performance_simulation(Codeword,SNR,H,"ALTN",NREALS,MAX)
    end
    if TABL
        @time FER_tabl, BER_tabl, Iters_tabl = 
            performance_simulation(Codeword,SNR,H,"TABL",NREALS,MAX)
    end
    if MSUM
        @time FER_msum, BER_msum, Iters_msum = 
            performance_simulation(Codeword,SNR,H,"MSUM",NREALS,MAX)
    end
    if  LBP
        @time FER_lbp, BER_lbp, Iters_lbp = 
            performance_simulation(Codeword,SNR,H, "LBP",NREALS,MAX)
    end
    if iLBP
        @time FER_ilbp, BER_ilbp, Iters_ilbp = 
            performance_simulation(Codeword,SNR,H,"iLBP",NREALS,MAX)
    end
    if  RBP
        @time FER_rbp, BER_rbp, Iters_rbp = 
            performance_simulation(Codeword,SNR,H,"RBP",NREALS,MAX_RBP)
    end
    if LRBP
        @time FER_lrbp, BER_lrbp, Iters_lrbp = 
            performance_simulation(Codeword,SNR,H,"LRBP",NREALS,MAX_RBP)
    end
    ################################### PLOTTING ###################################
    plotlyjs()
    lim = log10(1/NREALS)
    yaxis = Vector{Vector{Float64}}(undef,0)
    fer_labels = Vector{String}(undef,0)
    if MKAY
        append!(yaxis,[FER])
        push!(fer_labels,"SPA FL (McKay)")
    end
    if TANH
        append!(yaxis,[FER_tanh])
        push!(fer_labels,"SPA FL (tanh)")
    end
    if ALTN
        append!(yaxis,[FER_altn])
        push!(fer_labels,"SPA FL (alt)")
    end
    if TABL
        append!(yaxis,[FER_tabl])
        push!(fer_labels,"SPA FL (table)")
    end
    if MSUM
        append!(yaxis,[FER_msum])
        push!(fer_labels,"MIN SUM")
    end
    if LBP
        append!(yaxis,[FER_lbp])
        push!(fer_labels,"SPA LBP")
    end
    if iLBP
        append!(yaxis,[FER_ilbp])
        push!(fer_labels,"SPA iLBP")
    end
    if RBP
        append!(yaxis,[FER_rbp])
        push!(fer_labels,"SPA RBP")
    end
    if LRBP
        append!(yaxis,[FER_lrbp])
        push!(fer_labels,"SPA LRBP")
    end
    fer_labels = permutedims(fer_labels)

    ber_labels = Vector{String}(undef,0)
    for snr in SNR_db
        push!(ber_labels,"SNR (dB) = $snr")
    end
    ber_labels = permutedims(ber_labels)

    display(
        plot(
            SNR_db,yaxis,
            xlabel="SNR",
            label=fer_labels,
            lw=2,
            title="FER (girth = $girth)",
            ylims=(lim,0)
        )
    )

    if PLOT_BER
        if MKAY
            display(
                plot(
                    1:MAX,
                    BER,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (McKay)",
                    ylims=(lim-1,0)
                )
            )
        end
        if TANH
            display(
                plot(
                    1:MAX,
                    BER_tanh,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (tanh)",
                    ylims=(lim-1,0)
                )
            )
        end
        if ALTN
            display(
                plot(
                    1:MAX,
                    BER_altn,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (alt)",
                    ylims=(lim-1,0)
                )
            )
        end
        if TABL
            display(
                plot(
                    1:MAX,
                    BER_tabl,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (table)",
                    ylims=(lim-1,0)
                )
            )
        end
        if MSUM
            display(
                plot(
                    1:MAX,
                    BER_msum,
                    label=ber_labels,
                    lw=2,
                    title="BER MIN SUM",
                    ylims=(lim-1,0)
                )
            )
        end
        if  LBP
            display(
                plot(
                    1:MAX,
                    BER_lbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA LBP",
                    ylims=(lim-1,0)
                )
            )
        end
        if iLBP
            display(
                plot(
                    1:MAX,
                    BER_ilbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA iLBP",
                    ylims=(lim-1,0)
                )
            )
        end
        if  RBP
            display(
                plot(
                    1:MAX_RBP,
                    BER_rbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA RBP (decay cte = $DECAYCTE)",
                    ylims=(lim-1,0)
                )
            )
        end
        if LRBP
            display(
                plot(
                    1:MAX_RBP,
                    BER_lrbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA LRBP (decay cte = $DECAYCTE)",
                    ylims=(lim-1,0)
                )
            )
        end
    end

    if HISTOGRAMS
        for i in eachindex(SNR)
            display(
                histogram(
                    [Iters_tanh[i,:] Iters_ilbp[i,:]],
                    layout=grid(2,1),
                    xlims=(0,MAX+1),
                    labels=["Flooding" "LRBP"],
                    title="SNR (dB) = $(SNR_db[i])"
                )
            )
        end
    end
end