
# transform EbN0 in standard deviations
STDEV_PROF = sqrt.(exp10.(-EbN0_TEST/10) / (2*RR))

print_algorithm_details(1,ALGO_PROF,BPTYPE,DECAY_TEST)

@time @profview simcore(
        KK,
        CODE_LENGTH,
        G_CRC,
        STDEV_PROF,
        HH,
        PP,
        NC,
        NV,
        E_H,
        PROTOCOL,
        LIFTSIZE,
        ALGO_PROF,
        BPTYPE,
        MAX_FRAME_ERRORS,
        MAXITER_TEST,
        RAYL,
        C_DR_ITER,
        DECAY_TEST,
        LISTSIZES,
        CI_GAMMA,
        RGN_SEEDS[1],
        TEST,
        false,
        PAR_PAR
    )