################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# LDPC coding-decoding performance simulation using Believe Propagation (BP)

################################ JULIA PACKAGES ################################

using LinearAlgebra
using Random
using Statistics
using Plots
using SparseArrays
using CSV, DataFrames
using Dates
using DelimitedFiles

##################################### SAVE #####################################

if length(ARGS) == 0
    SAVE = false
elseif ARGS[1] == "true"
    SAVE = true
    NOW = string(now())
    mkdir("./Saved Data/"*NOW)
    FILE = open("./Saved Data/"*NOW*"/output.txt", "w")
else
    SAVE = false
end

################################## CONSTANTS ###################################

const INF = typemax(Int64)
const INFFLOAT = 1e3
const NINFFLOAT = -INFFLOAT
const ALPHA = 0.7              # Min-Sum attenuation factor
const TABLESIZE = 8192
const TABLERANGE = 10
const SIZE_PER_RANGE = TABLESIZE/TABLERANGE

# Random seed
SEED::Int = 1111

################################ CONTROL FLAGS #################################

TEST::Bool = true
PRIN::Bool = false
PROF::Bool = false
STOP::Bool = false # stop simulation at zero syndrome (if true, BER curves are
# not printed)

################################## PARAMETERS ##################################

MAXITER::Int = 50
# FACTORS = [0.7, 0.8, 0.9, 1.0]
# FACTORS = collect(0.1:0.1:1.0)
FACTORS = [0.8, 0.9, 1.0]
EbN0 = [1.25]
# EbN0 = [1.2]
TRIALS = [256000]
# TRIALS = [128, 1280, 12800, 128000]

# TEST
MAXITER_TEST::Int = 4
EbN0_TEST::Float64 = 1.5
TRIALS_TEST::Int = 1
DECAY_TEST::Float64 = 1.0

################################### SCHEDULE ###################################

MODES = ["Flooding","LBP","RBP","RBP relative","List-RBP","NW-RBP","VN-RBP","VN-RBP-ALT","List-VN-RBP"]
NUM_MODES = length(MODES)
ACTIVE = zeros(Bool,NUM_MODES)
LISTSIZES = zeros(Int,4)

# BP type: "MKAY", "TANH", "FAST", "TABL", "MSUM"
BPTYPES = Vector{String}(undef,NUM_MODES)

# maximum number of BP iterations
MAXITERS = zeros(Int,NUM_MODES)

DECAYS = Vector{Vector{<:AbstractFloat}}(undef,NUM_MODES)
for i in 1:NUM_MODES
    DECAYS[i] = [0.0]
end

i = 1
# Flooding
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# LBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER

# RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# RBP relative
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List-RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# NW-RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Variable Node RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# Variable Node RBP ALT
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List Variable Node RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = "FAST"
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# List sizes (min values = 4 and 2)
LISTSIZES[1] = 32
LISTSIZES[2] = 2

######################## CODE LENGTH, RATE AND PROTOCOL ########################

# Transmitted message length
GG::Int = 2064
# Effective Rate
RR::Float64 = 1/3                       # WiMAX compatibility offset
# LDPC protocol: NR5G = NR-LDPC (5G); PEG = PEG; WiMAX = IEEE80216e;
PROTOCOL::String = "NR5G"
    LAMBDA = [0.21, 0.25, 0.25, 0.29, 0]
    RO = [1.0, 0, 0, 0, 0, 0]

include("setup.jl")
if !TEST
    include("plot_or_save.jl")
end

### WiMAX: N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
                                     # 1632,1728,1824,1920,2016,2112,2208,2304}.    
                  # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.