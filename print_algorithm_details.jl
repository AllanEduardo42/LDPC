function print_algorithm_details(
    count::Int,
    algorithm::String,
    bptype::String,
    decay::Float64,
)

    str = """
________________________________________________________________________________

                                  Algorithm $count                                  
________________________________________________________________________________

Message passing algorithm: $algorithm (using """
    if bptype == "MKAY"
        str *="Mckay's SPA method)"
    elseif bptype == "TANH"
        str *="LLR-SPA calculated by tanh function)"
    elseif bptype == "TABL"
        str *="LLR-SPA precalculated in look-up table)"
    elseif bptype == "MSUM"
        str *="LLRs calculated by min-sum algorithm)"
    end
    if decay != 0.0
        str *= "\nResidual decaying factor: $decay"
    end
    if algorithm == "List-RBP"
        str *= "\nList 1 size: $(LISTSIZES[1])\nList 2 size: $(LISTSIZES[2])"
    elseif algorithm == "C&DR-RBP"
        str *= "\nReturn activation at iteration $(C_DR_ITER)"
    end
    str *= "\n"
    println(str)
    if SAVE
        println(FILE,str)
    end
end