################################################################################
# Allan Eduardo Feitosa
# 3 out 2024
# First layer of the routine to estimate the LPDC performance (FER BER x SNR)

include("simcore.jl")

function 
    performance_sim(
        snr::Vector{<:Real},
        mode::String,
        trials::Integer,
        maxiter::Integer;
        rgn_samples_seeds=ones(Int,NTHREADS),
        printtest=false    
    )

############################### CHECK VALID MODE ###############################
    if !(FLOOTYPE == "MKAY" || FLOOTYPE == "TANH" || FLOOTYPE == "ALTN" || 
        FLOOTYPE == "TABL" || FLOOTYPE == "MSUM")      
        throw(
            ArgumentError(
                "$FLOOTYPE is not a valid flooding mode"
            )
        )
    end    

    if mode == "Flooding"
        supermode, mode = mode, FLOOTYPE        
    elseif mode == "LBP" || mode == "iLBP"
        supermode = "LBP"    
    elseif mode == "RBP" || mode == "Local-RBP" || mode == "Random-RBP"||
           mode == "List-RBP" 
        supermode = "RBP"
    else
        throw(
            ArgumentError(
                "$mode is not a valid mode"
            )
        )
    end

########################### PRINT SIMULATION DETAILS ###########################
    
    # if trials = 1, set test mode
    test = (trials < 3) ? true : false

    if supermode == "RBP"
        if mode == "RBP"
            rbpfactor = DECAYRBP
        elseif mode == "Local-RBP"
            rbpfactor = DECAYLRBP
        elseif mode == "Random-RBP"
            rbpfactor = DECAYRRBP
        elseif mode == "List-RBP"
            rbpfactor = DECAYLIST
        end           
    else
        rbpfactor = nothing
    end    
    
    if test && printtest
        println()
        print("###################### Starting simulation (Testing mode) #####")
        println("#################")
        println()
    elseif !test
        println()
        print("############################# Starting simulation #############")
        println("#################")
        println()
        println("Number of trials: $trials")
    end
    if !test || printtest
        if supermode == "Flooding"
            print("Message passing protocol: Flooding (using ")
            if FLOOTYPE == "MKAY"
                println("Mckay's SPA method)")
            elseif FLOOTYPE == "TANH"
                println("LLR-SPA calculated by tanh)")
            elseif FLOOTYPE == "ALTN"
                println("LLR-SPA calculated by ϕ function)")
            elseif FLOOTYPE == "TABL"
                println("LLR-SPA precalculated in look-up table)")
            elseif FLOOTYPE == "MSUM"
                println("LLRs calculated by min-sum algorithm)")
            end
        elseif supermode == "RBP"
            print("Message passing protocol: RBP (residuals calculated by ")
            if FLOOTYPE == "MKAY"
                println("Mckay's SPA method)")
            elseif FLOOTYPE == "TANH"
                println("tanh)")
            elseif FLOOTYPE == "ALTN"
                println("ϕ function)")
            elseif FLOOTYPE == "TABL"
                println("look-up table)")
            elseif FLOOTYPE == "MSUM"
                println("min-sum algorithm)")
            end
        else
            println("Message passing protocol: $mode")
        end

        println("Maximum number of iterations: $maxiter")
        println("Simulated for SNR (dB): $snr")
        println("Stop at zero syndrome ? $STOP")
        (supermode == "RBP") ? println("Decaying factor: $rbpfactor") : nothing
        (mode == "Random-RBP") ? println("Sample size: $SAMPLESIZE") : nothing 
        println()
    end

################################ MULTITHREADING ################################

    trials_multh = trials ÷ NTHREADS

    if !test
        K = length(SNR)
        if MTHR
            decoded, ber = zeros(maxiter,K,NTHREADS), zeros(maxiter,K,NTHREADS)
        else
            decoded, ber = zeros(maxiter,K), zeros(maxiter,K)
        end
        for k in 1:K
            if MTHR
                Threads.@threads for i in 1:NTHREADS
                    decoded[:,k,i], ber[:,k,i] = performance_simcore(
                                            Msg,
                                            Cword,
                                            snr[k],
                                            H,
                                            Zc,
                                            mode,
                                            supermode,
                                            FLOOTYPE,
                                            ALPHA,
                                            trials_multh,
                                            maxiter,
                                            STOP,
                                            FAST,
                                            rbpfactor,
                                            LISTSIZE,
                                            LISTSIZE2,
                                            SAMPLESIZE,
                                            Rgn_noise_seeds[i],
                                            rgn_samples_seeds[i],
                                            test,
                                            printtest)
                end
            else
                decoded[:,k], ber[:,k] = performance_simcore(
                                            Msg,
                                            Cword,
                                            snr[k],
                                            H,
                                            Zc,
                                            mode,
                                            supermode,
                                            FLOOTYPE,
                                            ALPHA,
                                            trials,
                                            maxiter,
                                            STOP,
                                            FAST,
                                            rbpfactor,
                                            LISTSIZE,
                                            LISTSIZE2,
                                            SAMPLESIZE,
                                            Rgn_noise_seeds[1],
                                            rgn_samples_seeds[1],
                                            test,
                                            printtest)
            end
        end

        if MTHR
            FER = zeros(maxiter,K)
            BER = zeros(maxiter,K)
            FER .= 1 .- sum(decoded,dims=3)/trials
            BER .= sum(ber,dims=3)/(trials*N)
        else
            FER = 1 .- decoded/trials
            BER = ber/(trials*N)
        end

        return log10.(FER), log10.(BER)
    
    else # IF TESTING

        Lr, Lq = performance_simcore(
                                    Msg,
                                    Cword,
                                    snr[1],
                                    H,
                                    Zc,
                                    mode,
                                    supermode,
                                    FLOOTYPE,
                                    ALPHA,
                                    trials,
                                    maxiter,
                                    STOP,
                                    FAST,
                                    rbpfactor,
                                    LISTSIZE,
                                    LISTSIZE2,
                                    SAMPLESIZE,
                                    Rgn_noise_seeds[1],
                                    rgn_samples_seeds[1],
                                    test,
                                    # testsignal,
                                    printtest)
        
        return Lr, Lq        
    end
end