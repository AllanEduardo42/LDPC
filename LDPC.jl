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
#using CSV, DataFrames
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
const MAXLR = 1e3
const MINLR = -MAXLR
# const ALPHA = 0.73                     # Min-Sum attenuation factor
# const ALPHA = 0.78
const ALPHA = 0.755
# const ALPHA2 = 0.86 
# const ALPHA2 = 0.91
const ALPHA2 = 0.885
const TABLESIZE = 8192
const TABLERANGE = 10
const SIZE_PER_RANGE = TABLESIZE/TABLERANGE

# Random seed
SEED::Int = 1111

################################ CONTROL FLAGS #################################

TEST::Bool = true
PRIN::Bool = true
PROF::Bool = false
STOP::Bool = true # stop simulation at zero syndrome (if true, BER curves are
# not printed)

############################### TEST PARAMETERS ################################

# TEST
### Maximum number of BP iterations
MAXITER_TEST::Int = 1
### EbN0
EbN0_TEST::Float64 = 2.0
### Number of Monte Carlo Trials
TRIALS_TEST::Int = 1
### Residual Decay factors
DECAY_TEST::Float64 = 1.0

################################## PARAMETERS ##################################

### Maximum number of BP iterations
MAXITER::Int = 50

### EbN0
# EbN0 = [1.0]
EbN0 = [1.0, 1.5, 2.0, 2.5]
# EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]

### Number of Monte Carlo Trials
# TRIALS = [1280]
TRIALS = [128, 1280, 12800, 128000]
# TRIALS = [128, 1280, 12800, 128000, 1280000]

### Residual Decay factors
FACTORS = [1.0]
# FACTORS = [0.7, 0.8, 0.9, 1.0]

############################### LDPC ALGORITHMS ################################

ALGORITHMS = ["Flooding",        # Flooding
              "LBP",             # Layered Believe Propagation
              "RBP",             # Residual Believe Propagation
              "RD-RBP",          # Residual Decay RBP
              "NW-RBP",          # Node Wise RBP
              "SVNF",            # Silent Variable Node Free RBP
              "D-SVNF",          # Dynamic SVNF
              "List-RBP",        # List RBP
              "C-RBP",           # Consensus RBP
              "C&R-RBP",         # Consensus & Return RBP
              "C&DR-RBP",        # Consensus & Delayed Return RBP
              "VC-RBP",          # Variable to Check RBP
              "OV-RBP"           # Oscillating Variable Node RBP
              ]

NUM_MODES = length(ALGORITHMS)
ACTIVE = zeros(Bool,NUM_MODES)

# BP type: "MKAY", "TANH", "TABL", "MSUM", "MSUM2"
BPTYPES = Vector{Vector{String}}(undef,NUM_MODES)

# maximum number of BP iterations
MAXITERS = zeros(Int,NUM_MODES)

DECAYS = Vector{Vector{<:AbstractFloat}}(undef,NUM_MODES)
for i in 1:NUM_MODES
    DECAYS[i] = [0.0]
end

### Flag to activate all algorithms
ACTIVE_ALL = false

i = 1
# Flooding
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER

# LBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER

# RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER

# RD-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# NW-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER

# SVNF
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER

# D-SVNF
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER

# List-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS
# List sizes (min values = 4 and 2)
LISTSIZES = zeros(Int,2)
LISTSIZES[1] = 16
LISTSIZES[2] = 2

# C-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# C&R-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# C&DR-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# VC-RBP
i += 1
ACTIVE[i] = 0
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

# OV-RBP
i += 1
ACTIVE[i] = 1
BPTYPES[i] = ["TANH"]
MAXITERS[i] = MAXITER
DECAYS[i] = FACTORS

######################## CODE LENGTH, RATE AND PROTOCOL ########################

# Transmitted message length
GG::Int = 576
# Effective Rate
RR::Float64 = 1/2                       # WiMAX compatibility offset
# LDPC protocol: NR5G = NR-LDPC (5G); PEG = PEG; WiMAX = IEEE80216e;
PROTOCOL::String = "WiMAX"
    LAMBDA = [0.21, 0.25, 0.25, 0.29, 0]
    RO = [1.0, 0, 0, 0, 0, 0]

### WiMAX: N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
                                     # 1632,1728,1824,1920,2016,2112,2208,2304}.    
                  # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.

# NR5G: A takes values in
# [    24,   32,   40,   48,   56,   64,   72,   80,   88,   96,  104,  112,  
#     120,  128,  136,  144,  152,  160,  168,  176,  184,  192,  208,  224,  
#     240,  256,  272,  288,  304,  320,  336,  352,  368,  384,  408,  432,  
#     456,  480,  504,  528,  552,  576,  608,  640,  672,  704,  736,  768,  
#     808,  848,  888,  928,  984, 1032, 1064, 1128, 1160, 1192, 1224, 1256,  
#    1288, 1320, 1352, 1416, 1480, 1544, 1608, 1672, 1736, 1800, 1864, 1928,  
#    2024, 2088, 2152, 2216, 2280, 2408, 2472, 2536, 2600, 2664, 2728, 2792, 
#    2856, 2976, 3104, 3240, 3368, 3496, 3624, 3752, 3824]

#################################### SETUP #####################################
include("setup.jl")

################################# PLOT RESULTS #################################
if !TEST && !SAVE
    include("plot_results.jl")
end