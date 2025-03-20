################################################################################
# Allan Eduardo Feitosa
# 3 out 2024
# First layer of the routine to estimate the LPDC performance (FER BER x SNR)

include("simcore.jl")

function 
    performance_sim(
        snr::Union{Vector{<:AbstractFloat},AbstractFloat},
        mode::String,
        trials::Union{Vector{<:Integer},Integer},
        maxiter::Integer,
        bptype::String,
        decay::AbstractFloat 
    )

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
    println(str)
    if SAVE
        println(FILE,str)
    end
    str = "Number of trials: $trials"
    println(str)
    if SAVE
        println(FILE,str)
    end
    str = "Message passing protocol: $mode (using "
    if bptype == "MKAY"
        str = str*"Mckay's SPA method)"
    elseif bptype == "TANH"
        str = str*"LLR-SPA calculated by tanh)"
    elseif bptype == "FAST"
        str = str*"LLR-SPA calculated by fast tanh)"
    elseif bptype == "ALTN"
        str = str*"LLR-SPA calculated by ϕ function)"
    elseif bptype == "TABL"
        str = str*"LLR-SPA precalculated in look-up table)"
    elseif bptype == "MSUM"
        str = str*"LLRs calculated by min-sum algorithm)"
    end
    println(str)
    if SAVE
        println(FILE,str)
    end
    str=
    """Maximum number of iterations: $maxiter
    Number of threads (multithreading): $NTHREADS
    Simulated for SNR (dB): $snr
    Stop at zero syndrome ? $STOP"""
    println(str)
    if SAVE
        println(FILE,str)
    end
    if mode == "RBP" || mode == "Local-RBP" || mode == "List-RBP" ||
       mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        str = "RBP decaying factor: $decay"
        println(str)
        if SAVE
            println(FILE,str)
        end
    end
    if mode == "List-RBP" || mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        str = "List 1 size: $(LISTSIZES[1])\nList 2 size: $(LISTSIZES[2])"
        println(str)
        if SAVE
            println(FILE,str)
        end
    end
    println()
    if SAVE
        println(FILE,"")
    end

    ################################ MULTITHREADING ################################

    if TEST
        Lr, Lq, max_residues = simcore(
            AA,
            snr,
            HH,
            MM,
            NN,
            CN2VN,
            VN2CN,
            E_H,
            LDPC,
            ZF,
            NR_LDPC_DATA,
            mode,
            bptype,
            trials,
            maxiter,
            STOP,
            decay,
            LISTSIZES,
            THRES,
            RGN_NOISE_SEEDS[1],
            RGN_SAMPLE_SEEDS[1],
            RGN_MESSAGE_SEEDS[1];
            test=TEST,
            printtest = TEST ? PRIN : false)
        
        println()

        return Lr, Lq, max_residues

    else
        K = length(snr)
        sum_decoded = zeros(maxiter,K,NTHREADS)
        sum_ber = zeros(maxiter,K,NTHREADS)
        for k in 1:K
            stats = @timed Threads.@threads for i in 1:NTHREADS
                sum_decoded[:,k,i], sum_ber[:,k,i] = simcore(
                                                AA,
                                                snr[k],
                                                HH,
                                                MM,
                                                NN,
                                                CN2VN,
                                                VN2CN,
                                                E_H,
                                                LDPC,
                                                ZF,
                                                NR_LDPC_DATA,
                                                mode,
                                                bptype,
                                                trials[k]÷NTHREADS,
                                                maxiter,
                                                STOP,
                                                decay,
                                                LISTSIZES,
                                                THRES,
                                                RGN_NOISE_SEEDS[i],
                                                RGN_SAMPLE_SEEDS[i],
                                                RGN_MESSAGE_SEEDS[i])
            end
            str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
            println(str)
            if SAVE
                println(FILE,str)
            end
        end

        FER = sum(sum_decoded,dims=3)[:,:,1]
        for k = 1:K
            FER[:,k] ./= trials[k]
        end
        @. FER = 1 - FER
        BER = sum(sum_ber,dims=3)[:,:,1]
        for k = 1:K
            BER[:,k] ./= (NN*trials[k])
        end

        println()

        lowerfer = 1/maximum(trials)
        lowerber = lowerfer/NN
        replace!(x-> x < lowerfer ? lowerfer : x, FER)
        replace!(x-> x < lowerber ? lowerber : x, BER)

        return log10.(FER), log10.(BER)
        
    end

end