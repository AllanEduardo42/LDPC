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
MIN = true
LBP = true
RBP = true

PLOT_BER = true
HISTOGRAMS = false

################################ INCLUDED FILES ################################

# include("horizontal_update.jl")
include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")
include("auxiliary_functions.jl")
include("llr_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")
include("performance_estimation.jl")
include("test_SPA.jl")
include("lookupTable.jl")
include("SPA.jl")
include("PEG.jl")
include("GF2_functions.jl")
include("LBP.jl")
include("RBP.jl")

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1428
SEED2::Int64 = 5714

SIZE::Int64 = 1024
RANGE::Int64 = 20

SIZE_per_RANGE::Float64 = SIZE/RANGE

Phi = lookupTable()

NREALS::Int = 100
MAX::Int = 4

LR_idx = 9;

#################################### CODING ####################################

### PEG Algorithm

N = 512
M = 256
Random.seed!(SEED2)

D = rand([2,3],N)

@time H, girth = PEG!(D,M)

println("girth = ", girth)

@time G = gf2_nullspace(H)

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
    performance_estimation(
        C,
        Sigma,
        H,
        Checks2nodes,
        Nodes2checks,
        Phi,
        "TNH";
        nreals=1
    )
end
if ALT
    performance_estimation(
        C,
        Sigma,
        H,
        Checks2nodes,
        Nodes2checks,
        Phi,
        "ALT";
        nreals=1
    )
end
if TAB
    performance_estimation(
        C,
        Sigma,
        H,
        Checks2nodes,
        Nodes2checks,
        Phi,
        "TAB";
        nreals=1
    )
end
if MIN
    performance_estimation(
        C,
        Sigma,
        H,
        Checks2nodes,
        Nodes2checks,
        Phi,
        "MIN";
        nreals=1
    )
end
if LBP
    performance_estimation(
        C,
        Sigma,
        H,
        Checks2nodes,
        Nodes2checks,
        Phi,
        "LBP";
        nreals=1
    )
end
if RBP
    performance_estimation(
        C,
        Sigma,
        H,
        Checks2nodes,
        Nodes2checks,
        Phi,
        "RBP";
        nreals=1
    )
end
                             
########################### PERFORMANCE SIMULATION ############################

if TNH
    @time FER_tnh, BER_tnh, Iters_tnh, Lr_tnh = 
        performance_estimation(
            C,
            Sigma,
            H,
            Checks2nodes,
            Nodes2checks,
            Phi,
            "TNH";
            Lr_idx=LR_idx
        )
    ;
end
if ALT
    @time FER_alt, BER_alt, Iters_alt, Lr_alt = 
        performance_estimation(
            C,
            Sigma,
            H,
            Checks2nodes,
            Nodes2checks,
            Phi,
            "ALT";
            Lr_idx=LR_idx
        )
    ;
end
if TAB
    @time FER_tab, BER_tab, Iters_tab, Lr_tab = 
        performance_estimation(
            C,
            Sigma,
            H,
            Checks2nodes,
            Nodes2checks,
            Phi,
            "TAB";
            Lr_idx=LR_idx
        )
    ;
end
if MIN
    @time FER_min, BER_min, Iters_min, Lr_min = 
        performance_estimation(
            C,
            Sigma,
            H,
            Checks2nodes,
            Nodes2checks,
            Phi,
            "MIN";
            Lr_idx=LR_idx
        )
    ;
end
if LBP
    @time FER_lbp, BER_lbp, Iters_lbp, Lr_lbp = 
        performance_estimation(
            C,
            Sigma,
            H,
            Checks2nodes,
            Nodes2checks,
            Phi,
            "LBP";
            Lr_idx=LR_idx
        )
    ;
end
if RBP
    @time FER_rbp, BER_rbp, Iters_rbp, Lr_rbp = 
        performance_estimation(
            C,
            Sigma,
            H,
            Checks2nodes,
            Nodes2checks,
            Phi,
            "RBP";
            Lr_idx=LR_idx
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
                1:MAX,
                BER_rbp,
                label=ber_labels,
                lw=2,
                title="BER SPA RBP",
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