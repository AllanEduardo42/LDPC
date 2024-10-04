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
using CSV, DataFrames

################################ INCLUDED FILES ################################

include("performance_simulation.jl")
include("PEG.jl")
include("GF2_functions.jl")

################################ BP MODE FLAGS ################################

LBP::Bool = true
iLBP::Bool = true
 RBP::Bool = true
RRBP::Bool = true
LRBP::Bool = true

############################# FLOODING MODE FLAGS ##############################

MKAY::Bool = false
TANH::Bool = true
ALTN::Bool = false
TABL::Bool = false
MSUM::Bool = false
FAST::Bool = true # fast flooding update when using tanh mode (default:true)

################################ CONTROL FLAGS #################################

SAVE = true
PRINTING::Bool = false
PLOTBER::Bool = true
HISTOGRAMS::Bool = false

############################# SIMULATION CONSTANTS #############################

const INF = typemax(Int64)

# Seeds
SEED_NOISE::Int64 = 1428
SEED_GRAPH::Int64 = 5714
SEED_SAMPL::Int64 = 2857
SEED_MESSA::Int64 = 9999

# LookupTable
SIZE::Int64 = 1024
RANGE::Int64 = 20
SIZE_per_RANGE::Float64 = SIZE/RANGE

# Number of realizations and iterations
NREALS::Int = 100_000
MAX::Int = 30
MAXRBP::Int = 6
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are 
# not printed)

# decaying factor for RBP
DECAYRBP::Float64 = 0.9
DECAYLRBP::Float64 = 0.9
DECAYRRBP::Float64 = 0.9
SAMPLESIZE::Int = 51

##################################### SNR ######################################
SNRTEST = [4]
SNR = collect(1:1:4)

############################# PARITY-CHECK MATRIX #############################

### PEG Algorithm

N::Int64 = 512
M::Int64 = 256

rng_graph = Xoshiro(SEED_GRAPH)

D = rand(rng_graph,[2,3,4],N)

H, girth = PEG(D,M)

G = gf2_nullspace(H)

############################# MESSAGE AND CODEWORD #############################

rgn_message = Xoshiro(SEED_MESSA)

Message::Vector{Bool} = rand(rgn_message,Bool,N-M)

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

############################## JULIA COMPILATION ###############################
if MKAY
    Lr_mkay, Lq_mkay = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "MKAY",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if TANH
    Lr_tanh, Lq_tanh = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "TANH",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if ALTN
    Lr_altn, Lq_altn = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "ALTN",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if TABL
    Lr_tabl, Lq_tabl = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "TABL",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if MSUM
    Lr_msum, Lq_msum = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "MSUM",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if  LBP
      Lr_lbp, Lq_lbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "LBP",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if iLBP
    Lr_ilbp, Lq_ilbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "iLBP",
                                              1,
                                              MAX,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if RBP
    Lr_rbp, Lq_rbp =   performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "RBP",
                                              1,
                                              MAXRBP,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
if RRBP
    Lr_rrbp, Lq_rrbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "RRBP",
                                              1,
                                              MAXRBP,
                                              STOP,
                                              SEED_NOISE;
                                              rng_seed_sample=SEED_SAMPL,
                                              printing=PRINTING)
end
if LRBP
    Lr_lrbp, Lq_lrbp = performance_simulation(Codeword,
                                              SNRTEST,
                                              H,
                                              "LRBP",
                                              1,
                                              MAXRBP,
                                              STOP,
                                              SEED_NOISE;
                                              printing=PRINTING)
end
                             
############################ PERFORMANCE SIMULATION ############################
if NREALS > 1
    if MKAY
        @time FER_mkay, BER_mkay, Iters_mkay  = performance_simulation(
            Codeword,
            SNR,
            H,
            "MKAY",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if TANH
        @time FER_tanh, BER_tanh, Iters_tanh =performance_simulation(
            Codeword,
            SNR,
            H,
            "TANH",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if ALTN
        @time FER_altn, BER_altn, Iters_altn = performance_simulation(
            Codeword,
            SNR,
            H,
            "ALTN",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if TABL
        @time FER_tabl, BER_tabl, Iters_tabl = performance_simulation(
            Codeword,
            SNR,
            H,
            "TABL",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if MSUM
        @time FER_msum, BER_msum, Iters_msum = performance_simulation(
            Codeword,
            SNR,
            H,
            "MSUM",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if  LBP
        @time FER_lbp, BER_lbp, Iters_lbp = performance_simulation(
            Codeword,
            SNR,
            H,
            "LBP",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if iLBP
        @time FER_ilbp, BER_ilbp, Iters_ilbp = performance_simulation(
            Codeword,
            SNR,
            H,
            "iLBP",
            NREALS,
            MAX,
            STOP,
            SEED_NOISE)
    end
    if  RBP
        @time FER_rbp, BER_rbp, Iters_rbp = performance_simulation(
            Codeword,
            SNR,
            H,
            "RBP",
            NREALS,
            MAXRBP,
            STOP,
            SEED_NOISE)
    end
    if  RRBP
        @time FER_rrbp, BER_rrbp, Iters_rrbp = performance_simulation(
            Codeword,
            SNR,
            H,
            "RRBP",
            NREALS,
            MAXRBP,
            STOP,
            SEED_NOISE;
            rng_seed_sample=SEED_SAMPL)
    end
    if LRBP
        @time FER_lrbp, BER_lrbp, Iters_lrbp = performance_simulation(
            Codeword,
            SNR,
            H,
            "LRBP",
            NREALS,
            MAXRBP,
            STOP,
            SEED_NOISE)
    end
################################### PLOTTING ###################################
    plotlyjs()
    lim = log10(1/NREALS)
    yaxis = Vector{Vector{Float64}}(undef,0)
    fer_labels = Vector{String}(undef,0)
    if MKAY
        append!(yaxis,[FER_mkay])
        push!(fer_labels,"Flooding (McKay)")
    end
    if TANH
        append!(yaxis,[FER_tanh])
        push!(fer_labels,"Flooding (tanh)")
    end
    if ALTN
        append!(yaxis,[FER_altn])
        push!(fer_labels,"Flooding (alt)")
    end
    if TABL
        append!(yaxis,[FER_tabl])
        push!(fer_labels,"Flooding (table)")
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
        push!(fer_labels,"SPA RBP ($DECAYRBP)")
    end
    if RRBP
        append!(yaxis,[FER_rrbp])
        push!(fer_labels,"SPA RRBP ($DECAYRRBP)")
    end
    if LRBP
        append!(yaxis,[FER_lrbp])
        push!(fer_labels,"SPA LRBP ($DECAYLRBP)")
    end
    fer_labels = permutedims(fer_labels)

    ber_labels = Vector{String}(undef,0)
    for snr in SNR
        push!(ber_labels,"SNR (dB) = $snr")
    end
    ber_labels = permutedims(ber_labels)

    p = plot(
            SNR,yaxis,
            xlabel="SNR (dB)",
            label=fer_labels,
            lw=2,
            title="FER (Graph girth = $girth)",
            ylims=(lim,0)
        )

    display(p)
    SAVE ? savefig(p, "FER.png") : nothing

    if PLOTBER && !STOP
        if MKAY
            p = plot(
                1:MAX,
                BER_mkay,
                
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER Flooding (McKay)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_FLMKAY.png") : nothing
        end
        if TANH
            p = plot(
                1:MAX,
                BER_tanh,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER Flooding (tanh)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_FLTANH.png") : nothing
        end
        if ALTN
            p = plot(
                1:MAX,
                BER_altn,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER Flooding (alt)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_FLALTN.png") : nothing
        end
        if TABL
            p = plot(
                1:MAX,
                BER_tabl,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER Flooding (table)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_FLTABL.png") : nothing
        end
        if MSUM
            p = plot(
                1:MAX,
                BER_msum,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER Min-Sum",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_MINSUM.png") : nothing
        end
        if  LBP
            p = plot(
                1:MAX,
                BER_lbp,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER LBP",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_LBP.png") : nothing
        end
        if iLBP
            p = plot(
                1:MAX,
                BER_ilbp,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER iLBP",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_iLBP.png") : nothing
        end
        if  RBP
            p = plot(
                1:MAXRBP,
                BER_rbp,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER RBP (decay factor = $DECAYRBP)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_RBP.png") : nothing
        end
        if RRBP
            p = plot(
                1:MAXRBP,
                BER_rrbp,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER RRBP (decay factor = $DECAYRRBP)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_RRBP.png") : nothing
        end
        if LRBP
            p = plot(
                1:MAXRBP,
                BER_lrbp,
                xlabel="Iteration",
                label=ber_labels,
                lw=2,
                title="BER LRBP (decay factor = $DECAYLRBP)",
                ylims=(lim-2,0)
            )
            display(p)
            SAVE ? savefig(p,"BER_LRBP.png") : nothing
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
                    title="SNR (dB) = $(SNR)"
                )
            )
        end
    end

    if SAVE
    
        fers = []
        for i in eachindex(yaxis)
            push!(fers,(fer_labels[i],yaxis[i]))
        end
        FERS = Dict(fers)

        CSV.write("FERS.csv", DataFrame(FERS), header=true)

        bers = []
        padding = zeros(MAX-MAXRBP)
        if MKAY
            for i in eachindex(SNR)
                push!(bers,("Mckay (SNR=$(SNR))",BER_mkay[:,i]))
            end
        end
        if TANH
            for i in eachindex(SNR)
                push!(bers,("tanh (SNR=$(SNR))",BER_tanh[:,i]))
            end
        end
        if ALTN
            for i in eachindex(SNR)
                push!(bers,("altn (SNR=$(SNR))",BER_altn[:,i]))
            end
        end
        if TABL
            for i in eachindex(SNR)
                push!(bers,("tabl (SNR=$(SNR))",BER_tabl[:,i]))
            end
        end
        if MSUM
            for i in eachindex(SNR)
                push!(bers,("msum (SNR=$(SNR))",BER_msum[:,i]))
            end
        end
        if LBP
            for i in eachindex(SNR)
                push!(bers,("LBP (SNR=$(SNR))",BER_lbp[:,i]))
            end
        end
        if iLBP
            for i in eachindex(SNR)
                push!(bers,("iLBP (SNR=$(SNR))",BER_ilbp[:,i]))
            end
        end
        if RBP
            for i in eachindex(SNR)
                push!(bers,("RBP (SNR=$(SNR))",[BER_rbp[:,i];padding]))
            end
        end
        if RRBP
            for i in eachindex(SNR)
                push!(bers,("RRBP (SNR=$(SNR))",[BER_rrbp[:,i];padding]))
            end
        end
        if LRBP
            for i in eachindex(SNR)
                push!(bers,("RRBP (SNR=$(SNR))",[BER_lrbp[:,i];padding]))
            end
        end

        BERS = Dict(bers)

        CSV.write("BERS.csv", DataFrame(BERS), header=true)

    end

end