################################################################################
# Allan Eduardo Feitosa
# 3 out 2024
# First layer of the routine to estimate the LPDC performance (Fer Ber x SNR)

include("simcore.jl")

function 
    performance_sim(
        ebn0::Union{Vector{<:AbstractFloat},AbstractFloat},
        mode::String,
        trials::Union{Vector{<:Integer},Integer},
        maxiter::Integer,
        bptype::String,
        decay::AbstractFloat 
    )

    if bptype == "MKAY" && mode ≠ "Flooding"
        throw(error(lazy"The McKay method is only supported for Flooding"))
    elseif bptype == "TANH" && mode == "List-RBP"
        throw(error(lazy"The TANH method is not supported for List-RBP. Use the FAST, TABL or MSUM methods"))
    end

########################### PRINT SIMULATION DETAILS ###########################
    if TEST
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
        str *="LLR-SPA calculated by tanh)"
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
    Simulated for SNR (dB): $ebn0
    Stop at zero syndrome ? $STOP"""
    if mode == "RBP" || mode == "List-RBP" || mode == "VN-RBP" ||
       mode == "Genius-RBP" || mode == "NW-RBP"
        str *= "\nRBP decaying factor: $decay"
        str *= "\nRelative residues: $RELATIVE"
    end
    if mode == "List-RBP"
        str *= "\nList 1 size: $(LISTSIZES[1])\nList 2 size: $(LISTSIZES[2])"
    end
    str *= "\n"
    println(str)
    if SAVE
        println(FILE,str)
    end

    ################################ MULTITHREADING ################################

    if TEST
        @time Lr, Lq = simcore(
            AA,
            RR,
            GG,
            ebn0,
            HH,
            LL,
            UU,
            NC,
            NV,
            E_H,
            PROTOCOL,
            ZF,
            NR_LDPC_DATA,
            mode,
            bptype,
            trials,
            maxiter,
            STOP,
            decay,
            LISTSIZES,
            RELATIVE,
            RGN_NOISE_SEEDS[1],
            RGN_MESSAGE_SEEDS[1];
            test=TEST,
            printtest = TEST ? PRIN : false)
        
        println()
        if bptype == "TABL"
            Lr ./= SIZE_PER_RANGE
            Lq ./= SIZE_PER_RANGE
        end
        return Lr, Lq

    else
        K = length(ebn0)
        sum_decoded = zeros(Int,maxiter,K,NTHREADS)
        sum_ber = zeros(Int,maxiter,K,NTHREADS)
        for k in 1:K
            stats = @timed Threads.@threads for i in 1:NTHREADS
                sum_decoded[:,k,i], sum_ber[:,k,i] = simcore(
                                                AA,
                                                RR,
                                                GG,
                                                ebn0[k],
                                                HH,
                                                LL,
                                                UU,
                                                NC,
                                                NV,
                                                E_H,
                                                PROTOCOL,
                                                ZF,
                                                NR_LDPC_DATA,
                                                mode,
                                                bptype,
                                                trials[k]÷NTHREADS,
                                                maxiter,
                                                STOP,
                                                decay,
                                                LISTSIZES,
                                                RELATIVE,
                                                RGN_NOISE_SEEDS[i],
                                                RGN_MESSAGE_SEEDS[i])
            end
            str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
            println(str)
            if SAVE
                println(FILE,str)
            end
        end

        Fer = zeros(maxiter,K)
        Fer .= sum(sum_decoded,dims=3)
        for k = 1:K
            Fer[:,k] ./= trials[k]
        end
        @. Fer = 1 - Fer
        Ber = zeros(maxiter,K)
        Ber .= sum(sum_ber,dims=3)
        for k = 1:K
            Ber[:,k] ./= (NN*trials[k])
        end

        println()

        lowerfer = 1/maximum(trials)
        lowerber = lowerfer/NN
        replace!(x-> x < lowerfer ? lowerfer : x, Fer)
        replace!(x-> x < lowerber ? lowerber : x, Ber)

        # return log10.(Fer), log10.(Ber)
        return Fer, Ber

        return sum_decoded, sum_ber
        
    end

end