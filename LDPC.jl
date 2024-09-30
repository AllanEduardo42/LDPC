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

 LBP::Bool = true
iLBP::Bool = true
 RBP::Bool = false
LRBP::Bool = false

############################# FLOODING MODE FLAGS ##############################
MKAY::Bool = false
TANH::Bool = false
ALTN::Bool = false
TABL::Bool = false
MSUM::Bool = false

# fast flooding update when using tanh mode (default:true)
FAST::Bool = true

################################## TEST MODE ##################################
PRINTING::Bool = false
SNRTEST = [4]

####################### STOP WHEN SYNDROME IS ZERO FLAG ########################
STOP::Bool = false

################################### PLOTTING ###################################

PLOTBER::Bool = true
HISTOGRAMS::Bool = false

################################ INCLUDED FILES ################################

include("performance_simulation.jl")
include("PEG.jl")
include("GF2_functions.jl")

############################# SIMULATION CONSTANTS #############################

SEED_NOISE::Int64 = 1428
SEED_GRAPH::Int64 = 5714
SEED_SAMPL::Int64 = 2857
SEED_MESSA::Int64 = 9999

SIZE::Int64 = 1024
RANGE::Int64 = 20

SIZE_per_RANGE::Float64 = SIZE/RANGE

NREALS::Int = 1000
MAX::Int = 30
MAXRBP::Int = 5

DECAYRBP::Float64 = 0.9
DECAYLRBP::Float64 = 1
SAMPLESIZE::Int = 0

#################################### NOISE #####################################

SNR = collect(0:1:8)

############################# PARITY-CHECK MATRIX #############################

### PEG Algorithm

N::Int64 = 512
M::Int64 = 256

rng_graph = Xoshiro(SEED_GRAPH)

D = rand(rng_graph,[2,3,4],N)

H, girth = PEG!(D,M)

G = gf2_nullspace(H)

############################# MESSAGE AND CODEWORD #############################

rgn_message = Xoshiro(SEED_MESSA)

Message::Vector{Bool} = rand(rgn_message,Bool,N-M)

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
    Lr_mkay, Lq_mkay = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "MKAY",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING,
                                              stop=STOP)
end
if TANH
    Lr_tanh, Lq_tanh = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "TANH",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING, 
                                              stop=STOP)
end
if ALTN
    Lr_altn, Lq_altn = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "ALTN",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING, 
                                              stop=STOP)
end
if TABL
    Lr_tabl, Lq_tabl = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "TABL",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING,
                                              stop=STOP)
end
if MSUM
    Lr_msum, Lq_msum = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "MSUM",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING,
                                              stop=STOP)
end
if  LBP
      Lr_lbp, Lq_lbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "LBP",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING,
                                              stop=STOP)
end
if iLBP
    Lr_ilbp, Lq_ilbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "iLBP",
                                              1,
                                              MAX,
                                              SEED_NOISE;
                                              printing=PRINTING,
                                              stop=STOP)
end
if RBP
    Lr_rbp, Lq_rbp =   performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "RBP",
                                              1,
                                              MAXRBP,
                                              SEED_NOISE;
                                              rng_seed_sample=SEED_SAMPL,
                                              printing=PRINTING,
                                              stop=STOP)
end
if LRBP
    Lr_lrbp, Lq_lrbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "LRBP",
                                              1,
                                              MAXRBP,
                                              SEED_NOISE;
                                              printing=PRINTING,
                                              stop=STOP)
end
                             
############################ PERFORMANCE SIMULATION ############################
if NREALS > 1
    if MKAY
        @time FER_mkay, BER_mkay, Iters_mkay = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "MKAY",
                                   NREALS,
                                   MAX,
                                   SEED_NOISE;
                                   stop=STOP)
    end
    if TANH
        @time FER_tanh, BER_tanh, Iters_tanh = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "TANH",
                                   NREALS,
                                   MAX,
                                   SEED_NOISE;
                                   stop=STOP)
    end
    if ALTN
        @time FER_altn, BER_altn, Iters_altn = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "ALTN",
                                   NREALS,
                                   MAX,
                                   SEED_NOISE;
                                   stop=STOP)
    end
    if TABL
        @time FER_tabl, BER_tabl, Iters_tabl = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "TABL",
                                   NREALS,
                                   MAX,
                                   SEED_NOISE;
                                   stop=STOP)
    end
    if MSUM
        @time FER_msum, BER_msum, Iters_msum = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "MSUM",
                                   NREALS,
                                   MAX,
                                   SEED_NOISE;
                                   stop=STOP)
    end
    if  LBP
        @time FER_lbp, BER_lbp, Iters_lbp = 
            performance_simulation(Codeword,
                                  SNR,
                                  H,
                                  "LBP",
                                  NREALS,
                                  MAX,
                                  SEED_NOISE;
                                  stop=STOP)
    end
    if iLBP
        @time FER_ilbp, BER_ilbp, Iters_ilbp = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "iLBP",
                                   NREALS,
                                   MAX,
                                   SEED_NOISE;
                                   stop=STOP)
    end
    if  RBP
        @time FER_rbp, BER_rbp, Iters_rbp = 
            performance_simulation(Codeword,
                                   SNR,
                                   H,
                                   "RBP",
                                   NREALS,
                                   MAXRBP,
                                   SEED_NOISE;
                                   stop=STOP,
                                   rng_seed_sample=SEED_SAMPL)
    end
    if LRBP
        @time FER_lrbp, BER_lrbp, Iters_lrbp = 
            performance_simulation(Codeword,SNR,H,"LRBP",NREALS,MAXRBP,SEED_NOISE;
                                                                      stop=STOP)
    end
    ################################### PLOTTING ###################################
    plotlyjs()
    lim = log10(1/NREALS)
    yaxis = Vector{Vector{Float64}}(undef,0)
    fer_labels = Vector{String}(undef,0)
    if MKAY
        append!(yaxis,[FER_mkay])
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
    for snr in SNR
        push!(ber_labels,"SNR (dB) = $snr")
    end
    ber_labels = permutedims(ber_labels)

    display(
        plot(
            SNR,yaxis,
            xlabel="SNR (dB)",
            label=fer_labels,
            lw=2,
            title="FER (Graph girth = $girth)",
            ylims=(lim,0)
        )
    )

    if PLOTBER && !STOP
        if MKAY
            display(
                plot(
                    1:MAX,
                    BER_mkay,
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
                    1:MAXRBP,
                    BER_rbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA RBP (decay cte = $DECAYRBP)",
                    ylims=(lim-1,0)
                )
            )
        end
        if LRBP
            display(
                plot(
                    1:MAXRBP,
                    BER_lrbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA LRBP (decay cte = $DECAYLRBP)",
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
                    title="SNR (dB) = $(SNR[i])"
                )
            )
        end
    end
end