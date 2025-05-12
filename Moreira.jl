###############################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Test LDPC using Moreira's example

##################################### TEST #####################################

include("simcore.jl")
include("NR_LDPC_encode.jl")

const INF = typemax(Int64)
const INFFLOAT = 1e3
const NINFFLOAT = -INFFLOAT
const ALPHA = 0.875               # Min-Sum attenuation factor

HH = Matrix(Bool.(
     [0 1 0 1 0 1 1 1 0 0 0 1
      1 0 1 1 0 0 0 0 1 0 0 0
      0 1 0 0 1 0 1 0 0 0 0 1
      1 0 0 1 0 0 0 0 0 1 1 0
      0 0 1 0 1 1 0 0 0 1 0 0
      1 0 1 0 0 0 1 1 0 0 1 0
      0 1 0 0 0 1 0 1 1 1 0 0
      0 0 0 0 1 0 0 0 1 0 1 1]
      ))
    

MM,NN = size(HH)

NC  = make_cn2vn_list(HH)
NV  = make_vn2cn_list(HH)

GG = [gf2_mat_mult(gf2_inverse(HH[:,1:MM]),HH[:,MM+1:NN]); I(NN-MM)] 

### McKay Method

MSGTEST = [true, false, false, false]
CWORD = gf2_mat_mult(GG, MSGTEST) 
U = Float64.(2*CWORD .- 1) 

STDEV = 0.8
SNR = 10*log10(1/STDEV^2)

SIGNALTEST = [1.3129, 2.6584, 0.7413, 2.1745, 0.5981, −0.8323, −0.3962, −1.7586,
1.4905, 0.4084, −0.9290, 1.0765]

AA = length(MSGTEST)
RR = 1/2

PRINTTEST = true
MAXITER = 3
BPTYPE = "MKAY"

MLR, MLQ  = simcore(AA,
                    RR,
                    SNR,
                    HH,
                    GG,
                    NC,
                    NV,
                    nothing,
                    "PEG",
                    0,
                    nothing,
                    "Flooding",
                    BPTYPE,
                    1,
                    MAXITER,
                    false,
                    0.0,
                    [0],
                    false,
                    0,
                    0;
                    test=true,
                    printtest=PRINTTEST,
                    msgtest=MSGTEST,
                    noisetest=(SIGNALTEST-U)/STDEV)

# if BPTYPE == "MKAY"
#     LLR = log.(MLR[:,:,1]) - log.(MLR[:,:,2])
#     LLQ = log.(MLQ[:,:,1]) - log.(MLQ[:,:,2])
# else
#     LLR = MLR
#     LLQ = MLQ
# end
;

