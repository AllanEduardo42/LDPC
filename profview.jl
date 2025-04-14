################################################################################
# Allan Eduardo Feitosa
# 28 Feb 2025
# Script for profviewing


SNR_PROF = 1.2
Mode_prof = "RBP"
TRIALSPROF = 2^5
Maxiter_prof = 50
Bptypte_prof = "FAST"
Decay_prof = 0.9


@time @profview simcore(
    AA,
    SNR_PROF,
    HH,
    GG,
    CN2VN,
    VN2CN,
    E_H,
    PROTOCOL,
    ZF,
    NR_LDPC_DATA,
    Mode_prof,
    Bptypte_prof,
    TRIALSPROF,
    Maxiter_prof,
    STOP,
    Decay_prof,
    LISTSIZES,
    RELATIVE,
    RGN_NOISE_SEEDS[1],
    RGN_MESSAGE_SEEDS[1])