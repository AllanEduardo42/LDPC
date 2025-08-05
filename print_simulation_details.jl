function 
    print_simulation_details(
        test::Bool,
        trials::Vector{Int},
        mode::String,
        bptype::String,
        maxiter::Int,
        ebn0::Vector{Float64},
        decay::Float64,
    )
    if test
        str =
        """###################### Starting simulation (Testing mode) ######################
        """

    else
        str = 
        """############################# Starting simulation ##############################
        """
    end
    str *= """\nNumber of trials: $trials
    Message passing protocol: $mode (using """
    if bptype == "MKAY"
        str *="Mckay's SPA method)"
    elseif bptype == "TANH"
        str *="LLR-SPA calculated by tanh function)"
    elseif bptype == "TABL"
        str *="LLR-SPA precalculated in look-up table)"
    elseif bptype == "MSUM"
        str *="LLRs calculated by min-sum algorithm)"
    end
    str *=
    """\nMaximum number of iterations: $maxiter
    Number of threads (multithreading): $NTHREADS
    Eb/No (dB): $ebn0
    Stop at zero syndrome ? $STOP"""
    if mode == "RBP" || mode == "List-RBP" || mode == "List-VN-RBP" || mode == "VN-RBP" || mode == "TW-RBP" || mode == "C-RBP" || mode == "C-VN-RBP"
        str *= "\nRBP decaying factor: $decay"
        if mode == "List-RBP"
            str *= "\nList 1 size: $(LISTSIZES[1])\nList 2 size: $(LISTSIZES[2])"
        end
    end
    str *= "\n"
    println(str)
    if SAVE
        println(FILE,str)
    end
end