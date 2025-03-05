################################################################################
# Allan Eduardo Feitosa
# 28 Feb 2025
# Script for profviewing


SNR_PROF = 1.2
Mode_prof = "RBP"
TRIALSPROF = 2^5
Maxiter_prof = 2
Bptypte_prof = "FAST"
Decay_prof = 0.9


@profview simcore(
    A,
    SNR_PROF,
    H,
    E_H,
    LDPC,
    Zf,
    nr_ldpc_data,
    Mode_prof,
    Bptypte_prof,
    TRIALSPROF,
    Maxiter_prof,
    STOP,
    Decay_prof,
    Listsizes,
    Rgn_noise_seeds[1],
    Rgn_samples_seeds[1],
    Rgn_message_seeds[1])