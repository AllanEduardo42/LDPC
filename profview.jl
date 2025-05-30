################################################################################
# Allan Eduardo Feitosa
# 28 Feb 2025
# Script for profviewing

include("simcore.jl")

EbN0_PROF = 1.2
MODE_PROOF = "Flooding"
TRIALSPROF = 2^10
MAXITER_PROOF = 50
BPTYPE_PROOF = "FAST"
DECAY_PROOF = 0.9


@time @profview simcore(
    AA,
    RR,
    GG,
    EbN0_PROF,
    HH,
    LL,
    UU,
    NC,
    NV,
    E_H,
    PROTOCOL,
    ZF,
    NR_LDPC_DATA,
    MODE_PROOF,
    BPTYPE_PROOF,
    TRIALSPROF,
    MAXITER_PROOF,
    STOP,
    DECAY_PROOF,
    LISTSIZES,
    RELATIVE,
    RGN_NOISE_SEEDS[1],
    RGN_MESSAGE_SEEDS[1];
    test=false,
    printtest = false)