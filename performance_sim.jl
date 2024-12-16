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
        maxiter::Integer,
        bptype::String,
        decayfactor::AbstractFloat;
        printtest=false    
    )

########################### PRINT SIMULATION DETAILS ###########################
    
    # if trials = 1, set test mode
    test = (trials < 3) ? true : false  
    
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
        print("Message passing protocol: $mode (using ")
        if bptype == "MKAY"
            println("Mckay's SPA method)")
        elseif bptype == "TANH"
            println("LLR-SPA calculated by tanh)")
        elseif bptype == "FAST"
            println("LLR-SPA calculated by fast tanh)")
        elseif bptype == "ALTN"
            println("LLR-SPA calculated by ϕ function)")
        elseif bptype == "TABL"
            println("LLR-SPA precalculated in look-up table)")
        elseif bptype == "MSUM"
            println("LLRs calculated by min-sum algorithm)")
        end
        println("Maximum number of iterations: $maxiter")
        println("Simulated for SNR (dB): $snr")
        println("Stop at zero syndrome ? $STOP")
        println("Decaying factor: $decayfactor")
        (mode == "Random-RBP") ? println("Sample size: $SAMPLESIZE") : nothing
        (mode == "List-RBP") ? println("List 1 size: $LISTSIZE") : nothing 
        (mode == "List-RBP") ? println("List 2 size: $LISTSIZE2") : nothing 
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
                                            A,
                                            snr[k],
                                            H,
                                            E_H,
                                            zf,
                                            Zc,
                                            mode,
                                            bptype,
                                            trials_multh,
                                            maxiter,
                                            STOP,
                                            decayfactor,
                                            LISTSIZE,
                                            LISTSIZE2,
                                            Rgn_noise_seeds[i],
                                            Rgn_samples_seeds[i],
                                            Rgn_message_seeds[i],
                                            test,
                                            printtest)
                end
            else
                decoded[:,k], ber[:,k] = performance_simcore(
                                            A,
                                            snr[k],
                                            H,
                                            E_H,
                                            zf,
                                            Zc,
                                            mode,
                                            bptype,
                                            trials,
                                            maxiter,
                                            STOP,
                                            decayfactor,
                                            LISTSIZE,
                                            LISTSIZE2,
                                            Rgn_noise_seeds[1],
                                            Rgn_samples_seeds[1],
                                            Rgn_message_seeds[1],
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
                                    A,
                                    snr[1],
                                    H,
                                    E_H,
                                    zf,
                                    Zc,
                                    mode,
                                    bptype,
                                    trials,
                                    maxiter,
                                    STOP,
                                    decayfactor,
                                    LISTSIZE,
                                    LISTSIZE2,
                                    Rgn_noise_seeds[1],
                                    Rgn_samples_seeds[1],
                                    Rgn_message_seeds[1],
                                    test,
                                    printtest)
        
        return Lr, Lq        
    end
end