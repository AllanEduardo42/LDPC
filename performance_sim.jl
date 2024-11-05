################################################################################
# Allan Eduardo Feitosa
# 3 out 2024
# First layer of the routine to estimate the LPDC performance (FER BER x SNR)

include("performance_simcore.jl")

function 
    performance_sim(
        message::Vector{Bool},
        codeword::Vector{Bool},
        snr::Vector{<:Real},
        H::BitMatrix,
        mode::String,
        trials::Integer,
        maxiter::Integer,
        stop::Bool,
        rgn_noise_seeds::Vector{<:Integer};
        rgn_samples_seeds=ones(Int,NTHREADS),
        floomode="TANH",
        printtest=false    
    )

############################### CHECK VALID MODE ###############################
    if mode == "Flooding"
        if floomode == "MKAY" || floomode == "TANH" || floomode == "ALTN" || 
           floomode == "TABL" || floomode == "MSUM"      

            supermode, mode = mode, floomode
        else
            throw(
                ArgumentError(
                    "$floomode is not a valid flooding mode"
                )
            )
        end
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
            if floomode == "MKAY"
                println("Mckay's SPA method)")
            elseif floomode == "TANH"
                println("LLR-SPA calculated by tanh)")
            elseif floomode == "ALTN"
                println("LLR-SPA calculated by ϕ function)")
            elseif floomode == "TABL"
                println("LLR-SPA precalculated in look-up table)")
            elseif floomode == "MSUM"
                println("LLRs calculated by min-sum algorithm)")
            end
        else
            println("Message passing protocol: $mode")
        end

        println("Maximum number of iterations: $maxiter")
        println("Simulated for SNR (dB): $snr")
        println("Stop at zero syndrome ? $stop")
        (supermode == "RBP") ? println("Decaying factor: $rbpfactor") : nothing
        (mode == "Random-RBP") ? println("Sample size: $SAMPLESIZE") : nothing 
        println()
    end

################################ MULTITHREADING ################################

    trials_multh = trials ÷ NTHREADS

    if !test
        K = length(SNR)
        if MULTI
            decoded, ber = zeros(maxiter,K,NTHREADS), zeros(maxiter,K,NTHREADS)
        else
            decoded, ber = zeros(maxiter,K), zeros(maxiter,K)
        end
        for k in 1:K
            if MULTI
                Threads.@threads for i in 1:NTHREADS
                    decoded[:,k,i], ber[:,k,i] = performance_simcore(
                                            message,
                                            codeword,
                                            snr[k],
                                            H,
                                            mode,
                                            supermode,
                                            trials_multh,
                                            maxiter,
                                            stop,
                                            rbpfactor,
                                            rgn_noise_seeds[i],
                                            rgn_samples_seeds[i],
                                            test,
                                            printtest)
                end
            else
                decoded[:,k], ber[:,k] = performance_simcore(
                                            message,
                                            codeword,
                                            snr[k],
                                            H,
                                            mode,
                                            supermode,
                                            trials_multh,
                                            maxiter,
                                            stop,
                                            rbpfactor,
                                            rgn_noise_seeds[1],
                                            rgn_samples_seeds[1],
                                            test,
                                            printtest)
            end
        end

        if MULTI
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
                                    message,
                                    codeword,
                                    snr[1],
                                    H,
                                    mode,
                                    supermode,
                                    trials,
                                    maxiter,
                                    stop,
                                    rbpfactor,
                                    rgn_noise_seeds[1],
                                    rgn_samples_seeds[1],
                                    test,
                                    # testsignal,
                                    printtest)
        
        return Lr, Lq        
    end
end