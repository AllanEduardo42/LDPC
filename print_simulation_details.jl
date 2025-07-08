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
    elseif bptype == "RAW"
        str *="LLR-SPA calculated by raw-tanh)"
    elseif bptype == "FAST"
        str *="LLR-SPA calculated by fast tanh)"
    elseif bptype == "ALTN"
        str *="LLR-SPA calculated by ϕ function)"
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
    if mode == "RBP" || mode == "RBP relative" || mode == "List-RBP" || 
       mode == "VN-RBP" || mode == "Genius-RBP" || mode == "NW-RBP"
        str *= "\nRBP decaying factor: $decay"
    end
    if mode == "List-RBP"
        str *= "\nList 1 size: $(LISTSIZES[1])\nList 2 size: $(LISTSIZES[2])"
    end
    str *= "\n"
    println(str)
    if SAVE
        println(FILE,str)
    end
end