################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# LDPC coding-decoding performance simulation using Believe Propagation (BP)

############################## 1) JULIA PACKAGES ###############################

using LinearAlgebra
using Random
using Statistics
using Plots
using SparseArrays
using CSV, DataFrames

############################## 2) INCLUDED FILES ###############################

include("PEG.jl")
include("GF2_functions.jl")
include("IEEE80216e.jl")
include("NR_LDPC_encode.jl")
include("performance_sim.jl")
include("find_girth.jl")

################################# 3) CONSTANTS #################################

const INF = typemax(Int64)
const INFFLOAT = 1e2
const NINFFLOAT = -INFFLOAT

# Seeds
SEED_NOISE::Int = 1428
SEED_GRAPH::Int = 5714
SEED_SAMPL::Int = 2857
SEED_MESSA::Int = 9999

############################### 4) CONTROL FLAGS ###############################

MTHR::Bool = true                       
SAVE::Bool = false
PRIN::Bool = true
PLOT::Bool = true
STOP::Bool = true # stop simulation at zero syndrome (if true, BER curves are 
# not printed)

################################## 5) NUMBERS ##################################

TRIALS::Int = 10240*2
MAX::Int = 20
MAXRBP::Int = 10
SNRTEST = [4]
SNR = collect(1.0:0.4:2.2)

################################## 6) BP MODE ##################################

#Flooding
FLOO::Bool = true
    # Flooding type: "MKAY", "TANH", "ALTN", "TABL", "MSUM"
    FLOOTYPE = "TANH"
    # fast flooding update when using tanh mode (default:true)    
    FAST::Bool = true 

#LBP
_LBP::Bool = false 

#instantaneos-LBP
iLBP::Bool = false    

#RBP
_RBP::Bool = false 
    # RBP decay constant
    DECAYRBP::Float64 = 0.9  

#Random-RBP
RRBP::Bool = false  
    # Random-RBP decay constant  
    DECAYRRBP::Float64 = 0.5
    # Random-RBP sample size
    SAMPLESIZE::Int = 51
     
#Local-RBP
LRBP::Bool = false
    # Local-RBP decay constant    
    DECAYLRBP::Float64 = 0.5  

#List-RBP
LIST::Bool = true
    # List-RBP decay constant
    DECAYLIST::Float64 = 0.9
    # List-RBP size
    LISTSIZE::Int = fld(6384,128)
    LISTSIZE2::Int = 5

########################### 7) MESSAGE AND CODEWORD ############################

# Message (Payload) size
A::Int = 1008
# Rate
R::Float64 = 1/2
# LDPC protocol: 1 = NR-LDPC; 2 = PEG; 3 = IEEE80216e;
LDPC::Int = 3
    densities = 2:11

include("main.jl")