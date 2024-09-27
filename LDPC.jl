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

_LBP = false
ILBP = true
_RBP = false
LRBP = false

############################# FLOODING MODE FLAGS ##############################
FTNH = false
FALT = false
FTAB = false
FMSM = false

################################### PLOTTING ###################################

PLOT_BER = true
HISTOGRAMS = false

################################ INCLUDED FILES ################################

include("auxiliary_functions.jl")
include("performance_estimation.jl")
include("PEG.jl")
include("GF2_functions.jl")

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1428
SEED2::Int64 = 5714

SIZE::Int64 = 1024
RANGE::Int64 = 20

SIZE_per_RANGE::Float64 = SIZE/RANGE

NREALS::Int = 1000
MAX::Int = 30
MAX_RBP::Int = 5

LR_idx::Int = 9;
PENALTY::Float64 = 0.9

#################################### CODING ####################################

### PEG Algorithm

N::Int64 = 512
M::Int64 = 256
Random.seed!(SEED2)

D = rand([2,3,4],N)

H, girth = PEG!(D,M)

println("girth = ", girth)

G = gf2_nullspace(H)

K = size(G,2)

println("K = ", K)

Message = rand(Bool,K)

### codeword

C = gf2_mat_mult(Matrix(G), Message)

#################################### NOISE #####################################

SNR_db = collect(0:1:8)
SNR = exp10.(SNR_db/10)
Sigma = 1 ./ sqrt.(SNR)

############################# AUXILIARY CONSTANTS ##############################

Nodes2checks  = find_nodes2checks(H)
Checks2nodes  = find_checks2nodes(H)

############################## JULIA COMPILATION ###############################

if FTNH
    R, Lr_ftnh, Q, Lq_ftnh = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "FTNH",
        1,
        1;
        printing=false
    )
end
if FALT
    R, Lr_falt, Q, Lq_falt = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "FALT",
        1,
        1
    )
end
if FTAB
    R, Lr_ftab, Q, Lq_ftab = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "FTAB",
        1,
        1
    )
end
if FMSM
    R, Lr_msum, Q, Lq_msum = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "FMSM",
        1,
        1
    )
end
if _LBP
    R, Lr_lbp, Q, Lq_lbp = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "_LBP",
        1,
        1
    )
end
if _LBP
    R, Lr_ilbp, Q, Lq_ilbp = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "ILBP",
        1,
        1
    )
end
if _RBP
    R, Lr_rbp, Q, Lq_rbp, Edges = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "_RBP",
        1,
        1
    )
end
if LRBP
    R, Lr_lrbp, Q, Lq_lrbp, lEdges = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "LRBP",
        1,
        1
    )
end
                             
########################### PERFORMANCE SIMULATION ############################
if NREALS > 1
    if FTNH
        @time FER_ftnh, BER_ftnh, Iters_ftnh = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "FTNH",
                NREALS,
                MAX
            )
        ;
    end
    if FALT
        @time FER_falt, BER_falt, Iters_falt = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "FALT",
                NREALS,
                MAX
            )
        ;
    end
    if FTAB
        @time FER_ftab, BER_ftab, Iters_ftab = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "FTAB",
                NREALS,
                MAX
            )
        ;
    end
    if FMSM
        @time FER_msum, BER_msum, Iters_msum = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "FMSM",
                NREALS,
                MAX
            )
        ;
    end
    if _LBP
        @time FER_lbp, BER_lbp, Iters_lbp = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "_LBP",
                NREALS,
                MAX
            )
        ;
    end
    if ILBP
        @time FER_ilbp, BER_ilbp, Iters_ilbp = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "ILBP",
                NREALS,
                MAX
            )
        ;
    end
    if _RBP
        @time FER_rbp, BER_rbp, Iters_rbp = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "_RBP",
                NREALS,
                MAX_RBP
            )
        ;
    end
    if LRBP
        @time FER_lrbp, BER_lrbp, Iters_lrbp = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "LRBP",
                NREALS,
                MAX_RBP
            )
        ;
    end
    ################################### PLOTTING ###################################
    plotlyjs()
    lim = log10(1/NREALS)
    yaxis = Vector{Vector{Float64}}(undef,0)
    fer_labels = Vector{String}(undef,0)
    if FTNH
        append!(yaxis,[FER_ftnh])
        push!(fer_labels,"SPA FL (THN)")
    end
    if FALT
        append!(yaxis,[FER_falt])
        push!(fer_labels,"SPA FL (ALT)")
    end
    if FTAB
        append!(yaxis,[FER_ftab])
        push!(fer_labels,"SPA FL (TABLE)")
    end
    if FMSM
        append!(yaxis,[FER_msum])
        push!(fer_labels,"MIN SUM")
    end
    if _LBP
        append!(yaxis,[FER_lbp])
        push!(fer_labels,"SPA LBP")
    end
    if ILBP
        append!(yaxis,[FER_ilbp])
        push!(fer_labels,"SPA ILBP")
    end
    if _RBP
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
        if FTNH
            display(
                plot(
                    1:MAX,
                    BER_ftnh,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (TANH)",
                    ylims=(lim-1,0)
                )
            )
        end
        if FALT
            display(
                plot(
                    1:MAX,
                    BER_falt,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (ALT)",
                    ylims=(lim-1,0)
                )
            )
        end
        if FTAB
            display(
                plot(
                    1:MAX,
                    BER_ftab,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA FL (TABLE)",
                    ylims=(lim-1,0)
                )
            )
        end
        if FMSM
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
        if _LBP
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
        if ILBP
            display(
                plot(
                    1:MAX,
                    BER_ilbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA ILBP",
                    ylims=(lim-1,0)
                )
            )
        end
        if _RBP
            display(
                plot(
                    1:MAX_RBP,
                    BER_rbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA RBP (penalty = $PENALTY)",
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
                    title="BER SPA LRBP (penalty = $PENALTY)",
                    ylims=(lim-1,0)
                )
            )
        end
    end

    if HISTOGRAMS
        for i in eachindex(SNR)
            display(
                histogram(
                    [Iters_ftnh[i,:] Iters_ilbp[i,:]],
                    layout=grid(2,1),
                    xlims=(0,MAX+1),
                    labels=["Flooding" "local-RBP"],
                    title="SNR (dB) = $(SNR_db[i])"
                )
            )
        end
    end
end