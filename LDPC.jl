################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# LDPC coding-decoding performance simulation using Believe Propagation (BP)

################################ JULIA PACKAGES ################################

using LinearAlgebra
using Random
using LinearAlgebra
using Polynomials
using Statistics
using Plots
using SparseArrays
using Dates
using DelimitedFiles

################################## CONSTANTS ###################################

const INF = typemax(Int64)              
const MAXC2V = 1e3                      # saturate values for C2V (Inf approx)
const MINC2V = -MAXC2V                  # -Inf approx
const ALPHA = 0.75                      # Min-Sum attenuation factor #1
const ALPHA2 = 0.885                    # Min-Sum attenuation factor #2

# Random seed
SEED::Int = 1111

################################ CONTROL FLAGS #################################

TEST::Bool = true                       # Testing mode (few trials)
PRIN::Bool = true                       # Print info is in testing mode
PROF::Bool = false                      # profview
RAYL::Bool = false                      # Rayleigh fading channel

############################### TEST PARAMETERS ################################

### Maximum number of BP iterations
MAXITER_TEST::Int = 1
### EbN0
EbN0_TEST::Float64 = 2.0
### Maximum number of Frame Errors
ERRORS_TEST::Int = 1
### Residual Decay factors
DECAY_TEST::Float64 = 1.0
### CI-RBP gamma constant
CI_GAMMA = 0.15

################################## PARAMETERS ##################################

### Maximum number of Belief Propagation (BP) iterations
MAXITER::Int = 20

### EbN0
# EbN0 = [1.0, 1.25, 1.5, 1.75, 2.0]
EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
# EbN0 = [1.0]

### Maximum number of Frame Errors (at the last iteration)
ERRORS = 36*7

### Residual Decay factors
FACTORS = [0.85]
# FACTORS = [0.7, 0.8, 0.85, 0.9, 1.0]

############################### LDPC ALGORITHMS ################################

ALGORITHMS = ["Flooding",        # Flooding
              "LBP",             # Layered Believe Propagation
              "RBP",             # Residual Believe Propagation
              "RD-RBP",          # Residual Decay RBP
              "NW-RBP",          # Node Wise RBP
              "SVNF",            # Silent Variable Node Free RBP
              "List-RBP",        # List-RBP
              "C-RBP",           # Consensus RBP
              "C&R-RBP",         # Consensus & Return RBP
              "C&DR-RBP",        # Consensus & Delayed Return RBP
              "VC-RBP",          # Variable to Check RBP
              "OV-RBP",          # Oscillating Variable Node RBP
              "CI-RBP",          # Conditional Innovation RBP
              "UBP-RBP",         # Update Before Propagate RBP
              "RBP-D1VN"         # RBP with D1VN Processing Scheme
              ]
NUM_MODES = length(ALGORITHMS)

# Vector that indicates if each algorithm is active for simulation
# if ACTIVE[ALGO] == true, the performance of Algorithm[ALGO] will be simulated
ACTIVE = zeros(Bool,NUM_MODES)

# BP types: 
# "TANH": used the log-domain implementation of SPA (tanh functions)
# "MSUM": approximates the C2V functions using the min-sum algorithm
# "MSUMRBP": For RBP based algorithms, calculate only the residuals using min-sum
BPTYPES = Vector{Vector{String}}(undef,NUM_MODES)

# maximum number of BP iterations
MAXITERS = zeros(Int,NUM_MODES)

DECAYS = Vector{Vector{<:AbstractFloat}}(undef,NUM_MODES)
for ALGO in 1:NUM_MODES
    DECAYS[ALGO] = [0.0]
end

### Flag to activate all algorithms
ACTIVE_ALL = false

# For each algorithm:
# ACTIVE[ALGO] = 1 : the performance will be simulated
# BPTYPES[ALGO] = ["TANH", "MSUM"] : the BP types the algorithm will be simulated
# MAXITERS[ALGO] = 50 : number of BP iterations, defined for each active algorithm

ALGO = 1
# Flooding
ACTIVE[ALGO] = 1                           
BPTYPES[ALGO] = ["TANH"]                  
MAXITERS[ALGO] = MAXITER

# LBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# RD-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER
DECAYS[ALGO] = FACTORS

# NW-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# SVNF
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# List-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER
DECAYS[ALGO] = FACTORS
# List sizes (min values = 4 and 2)
LISTSIZES = zeros(Int,2)
LISTSIZES[1] = 16
LISTSIZES[2] = 2

# C-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER
DECAYS[ALGO] = FACTORS

# C&R-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER
DECAYS[ALGO] = FACTORS

# C&DR-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER
DECAYS[ALGO] = FACTORS
C_DR_ITER::Int = 4                           # Activation of Return in C&DR-RBP

# VC-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# OV-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# CI-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# UBP-RBP
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

# RBP-D1VN
ALGO += 1
ACTIVE[ALGO] = 0
BPTYPES[ALGO] = ["TANH"]
MAXITERS[ALGO] = MAXITER

######################## CODE LENGTH, RATE AND PROTOCOL ########################

# Transmitted message length
CODE_LENGTH::Int = 576
# Code Rate = RATE[1]/RATE[2]
RATE = [1, 2]              
# LDPC protocol: 5GNR = NR-LDPC (5G); PEG = PEG; WiMAX = IEEE80216e;
PROTOCOL::String = "5GNR"
    LAMBDA = [0.21, 0.25, 0.25, 0.29, 0]
    RO = [1.0, 0, 0, 0, 0, 0]

### WiMAX: N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
                                     # 1632,1728,1824,1920,2016,2112,2208,2304}.    
                  # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.

# 5GNR: A takes values in
# [    24,   32,   40,   48,   56,   64,   72,   80,   88,   96,  104,  112,  
#     120,  128,  136,  144,  152,  160,  168,  176,  184,  192,  208,  224,  
#     240,  256,  272,  288,  304,  320,  336,  352,  368,  384,  408,  432,  
#     456,  480,  504,  528,  552,  576,  608,  640,  672,  704,  736,  768,  
#     808,  848,  888,  928,  984, 1032, 1064, 1128, 1160, 1192, 1224, 1256,  
#    1288, 1320, 1352, 1416, 1480, 1544, 1608, 1672, 1736, 1800, 1864, 1928,  
#    2024, 2088, 2152, 2216, 2280, 2408, 2472, 2536, 2600, 2664, 2728, 2792, 
#    2856, 2976, 3104, 3240, 3368, 3496, 3624, 3752, 3824]

##################################### SAVE #####################################

if length(ARGS) == 0
    SAVE = false
elseif ARGS[1] == "true"
    SAVE = true
    NOW = string(now())
    DIRECTORY = "./Saved Data/"*NOW[1:10]*" "*NOW[12:16]*" $PROTOCOL $CODE_LENGTH $(RATE[1])|$(RATE[2])"
    for ebn0 in EbN0
        global DIRECTORY *= " $(ebn0)dB"
    end
    mkdir(DIRECTORY)
    FILE = open(DIRECTORY*"/output.txt", "w")
else
    SAVE = false
end

#################################### SETUP #####################################
include("setup.jl")

################################# PLOT RESULTS #################################
if !TEST && !SAVE
    include("plot_results.jl")
end

if SAVE
    println(FILE, "The End!")
    close(FILE)
end