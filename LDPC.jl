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

################################ SPA MODE FLAGS ################################

TNH = true
ALT = false
TAB = false
MIN = false
LBP = false
RBP = true
RBP_R = false

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

LR_idx = 9;

#################################### CODING ####################################

### PEG Algorithm

N::Int64 = 512
M::Int64 = 256
Random.seed!(SEED2)

D = rand([2,3,4],N)

H, girth = PEG!(D,M)

N_MESSAGES::Int64 = sum(H)

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

if TNH
    R_tnh, Lr_tnh, Q_tnh, LQ_tnh = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "TNH",
        1,
        MAX
    )
end
if ALT
    R_alt, Lr_alt, Q_alt, LQ_alt = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "ALT",
        1,
        MAX
    )
end
if TAB
    R_tab, Lr_tab, Q_tab, LQ_tab = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "TAB",
        1,
        MAX
    )
end
if MIN
    R_min, Lr_min, Q_min, LQ_min = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "MIN",
        1,
        MAX
    )
end
if LBP
    R_lbp, Lr_lbp, Q_lbp, LQ_lbp = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "LBP",
        1,
        MAX
    )
end
if RBP
    R_rbp, Lr_rbp, Q_rbp, LQ_rbp = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "RBP",
        1,
        MAX_RBP
    )
end
if RBP_R
    R_rbpr, Lr_rbpr, Q_rbpr, LQ_rbpr = performance_estimation(
        C,
        [Sigma[LR_idx]],
        H,
        Checks2nodes,
        Nodes2checks,
        "RBP_R",
        1,
        MAX_RBP
    )
end
                             
########################### PERFORMANCE SIMULATION ############################
if NREALS > 1
    if TNH
        @time FER_tnh, BER_tnh, Iters_tnh = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "TNH",
                NREALS,
                MAX
            )
        ;
    end
    if ALT
        @time FER_alt, BER_alt, Iters_alt = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "ALT",
                NREALS,
                MAX
            )
        ;
    end
    if TAB
        @time FER_tab, BER_tab, Iters_tab = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "TAB",
                NREALS,
                MAX
            )
        ;
    end
    if MIN
        @time FER_min, BER_min, Iters_min = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "MIN",
                NREALS,
                MAX
            )
        ;
    end
    if LBP
        @time FER_lbp, BER_lbp, Iters_lbp = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "LBP",
                NREALS,
                MAX
            )
        ;
    end
    if RBP
        @time FER_rbp, BER_rbp, Iters_rbp = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "RBP",
                NREALS,
                MAX_RBP
            )
        ;
    end
    if RBP_R
        @time FER_rbpr, BER_rbpr, Iters_rbpr = 
            performance_estimation(
                C,
                Sigma,
                H,
                Checks2nodes,
                Nodes2checks,
                "RBP_R",
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
    if TNH
        append!(yaxis,[FER_tnh])
        push!(fer_labels,"SPA TNH")
    end
    if ALT
        append!(yaxis,[FER_alt])
        push!(fer_labels,"SPA ALT")
    end
    if TAB
        append!(yaxis,[FER_tab])
        push!(fer_labels,"SPA TAB")
    end
    if MIN
        append!(yaxis,[FER_min])
        push!(fer_labels,"MIN SUM")
    end
    if LBP
        append!(yaxis,[FER_lbp])
        push!(fer_labels,"SPA LBP")
    end
    if RBP
        append!(yaxis,[FER_rbp])
        push!(fer_labels,"SPA RBP")
    end
    if RBP
        append!(yaxis,[FER_rbp])
        push!(fer_labels,"SPA RBP_R")
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
        if TNH
            display(
                plot(
                    1:MAX,
                    BER_tnh,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA TNH",
                    ylims=(lim-1,0)
                )
            )
        end
        if ALT
            display(
                plot(
                    1:MAX,
                    BER_alt,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA ALT",
                    ylims=(lim-1,0)
                )
            )
        end
        if TAB
            display(
                plot(
                    1:MAX,
                    BER_tab,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA TAB",
                    ylims=(lim-1,0)
                )
            )
        end
        if MIN
            display(
                plot(
                    1:MAX,
                    BER_min,
                    label=ber_labels,
                    lw=2,
                    title="BER MIN SUM",
                    ylims=(lim-1,0)
                )
            )
        end
        if LBP
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
        if RBP
            display(
                plot(
                    1:MAX_RBP,
                    BER_rbp,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA RBP",
                    ylims=(lim-1,0)
                )
            )
        end
        if RBP_R
            display(
                plot(
                    1:MAX_RBP,
                    BER_rbpr,
                    label=ber_labels,
                    lw=2,
                    title="BER SPA RBP_R",
                    ylims=(lim-1,0)
                )
            )
        end
    end

    if HISTOGRAMS
        for i in eachindex(SNR)
            display(
                histogram(
                    [Iters_alt[i,:] Iters_min[i,:]],
                    layout=grid(2,1),
                    xlims=(0,MAX+1),
                    labels=["SPA" "MIN SUM"],
                    title="SNR (dB) = $(SNR_db[i])"
                )
            )
        end
    end
end