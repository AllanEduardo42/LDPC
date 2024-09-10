################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# LDPC coding-decoding performance simulation

################################ JULIA PACKAGES ################################

using LinearAlgebra
using Random
using Statistics
using Plots

################################ SPA MODE FLAGS ################################

TNH = false
ALT = true
TAB = false
MIN = true

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

############################# SIMULATION CONSTANTS #############################

SEED::Int64 = 1427

SIZE::Int64 = 2048
RANGE::Int64 = 10

SIZE_per_RANGE::Float64 = SIZE/RANGE

Phi = lookupTable()

NREALS::Int = 1_000_0
MAX::Int = 10

#################################### CODING ####################################

### PEG Algorithm

N = 256
M = 128
D = rand([4,4],N)

@time H, girth = PEG!(D,M)

println("girth = ", girth)

@time G = GF2_nullspace(H)

K = size(G,2)

println("K = ", K)

Message = rand(Bool,K)

### codeword

C = bitwise_mat_mult(G, Message)

#################################### NOISE #####################################

SNR = collect(3:-0.3:0.3)
Sigma = 1 ./ SNR

############################# AUXILIARY CONSTANTS ##############################

Indices_col  = find_indices_col(H)
Indices_row  = find_indices_row(H)

############################## JULIA COMPILATION ###############################

if TNH
    performance_estimation(
        C,
        Sigma,
        H,
        Indices_row,
        Indices_col,
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
        Indices_row,
        Indices_col,
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
        Indices_row,
        Indices_col,
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
        Indices_row,
        Indices_col,
        Phi,
        "APP";
        nreals=1
    )
end
                             
########################### PERFORMANCE SIMULATION ############################

if TNH
    @time FER_tnh, BER_tnh, Iters_tnh = 
        performance_estimation(
            C,
            Sigma,
            H,
            Indices_row,
            Indices_col,
            Phi,
            "TNH"
        )
    ;
end
if ALT
    @time FER_alt, BER_alt, Iters_alt = 
        performance_estimation(
            C,
            Sigma,
            H,
            Indices_row,
            Indices_col,
            Phi,
            "ALT"
        )
    ;
end
if TAB
    @time FER_tab, BER_tab, Iters_tab = 
        performance_estimation(
            C,
            Sigma,
            H,
            Indices_row,
            Indices_col,
            Phi,
            "TAB"
        )
    ;
end
if MIN
    @time FER_app, BER_app, Iters_app = 
        performance_estimation(
            C,
            Sigma,
            H,
            Indices_row,
            Indices_col,
            Phi,
            "APP"
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
    append!(yaxis,[FER_app])
    push!(fer_labels,"MIN SUM")
end
fer_labels = permutedims(fer_labels)

ber_labels = Vector{String}(undef,0)
for snr in SNR
    push!(ber_labels,"SNR = $snr")
end
ber_labels = permutedims(ber_labels)

display(
    plot(
        SNR,yaxis,
        xlabel="SNR",
        label=fer_labels,
        lw=2,
        title="FER",
        ylims=(lim,0)
    )
)

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
    # plot!(girth*[1, 1]/2,[lim-1, 0],label="girth",lw=2,ls=:dot,lc=:black)
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
    # plot!(girth*[1, 1]/2,[lim-1, 0],label="girth",lw=2,ls=:dot,lc=:black)
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
    # plot!(girth*[1, 1]/2,[lim-1, 0],label="girth",lw=2,ls=:dot,lc=:black)
end
if MIN
    display(
        plot(
            1:MAX,
            BER_app,
            label=ber_labels,
            lw=2,
            title="BER SPA APP",
            ylims=(lim-1,0)
        )
    )
    # plot!(girth*[1, 1]/2,[lim-1, 0],label="girth",lw=2,ls=:dot,lc=:black)
end
