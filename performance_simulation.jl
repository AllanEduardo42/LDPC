################################################################################
# Allan Eduardo Feitosa
# 3 out 2024
# First layer of the routine to estimate the LPDC performance (FER BER x SNR)

include("performance_simulation_core.jl")

function 
    performance_simulation(
        codeword::Vector{Bool},
        snr::Vector{<:Real},
        H::BitMatrix,
        mode::String,
        nreals::Integer,
        max::Integer,
        stop::Bool,
        rgn_seed_noise::Integer;
        rng_seed_sample=1234,
        t_test=nothing,
        printing=false    
    )

############################### CHECK VALID MODE ###############################
    if mode == "MKAY" || mode == "TANH" || mode == "ALTN" || mode == "TABL" ||
        mode == "MSUM"
        
        supermode = "FLOO"

    elseif mode == "LBP" || mode == "iLBP"

        supermode = "LBP"
    
    elseif mode == "RBP" || mode == "LRBP" || mode == "RRBP"
        
        supermode = "RBP"

    else

        throw(
            ArgumentError(
                "$mode is not a valid mode"
            )
        )

    end

########################## PRINT SIMULATION DETAILS ############################
    
    # if nreals = 1, set test mode
    test = (nreals < 2) ? true : false
    if supermode == "RBP"
        if mode == "RBP"
            rbpfactor = DECAYRBP
        elseif mode == "LRBP"
            rbpfactor = DECAYLRBP
        elseif mode == "RRBP"
            rbpfactor = DECAYRRBP
        end           
    else
        rbpfactor = nothing
    end    
    
    if test && printing
        println()
        print("###################### Starting simulation (Testing mode) #####")
        println("#################")
        println()
    elseif !test
        println()
        print("############################# Starting simulation #############")
        println("#################")
        println()
        println("Number of trials: $nreals")
    end
    if !test || printing
        if supermode == "FLOO"
            print("Message passing protocol: Flooding (using ")
            if mode == "MKAY"
                println("Mckay's SPA method)")
            elseif mode == "TANH"
                println("LLR-SPA calculated by tanh)")
            elseif mode == "ALTN"
                println("LLR-SPA calculated by ϕ function)")
            elseif mode == "TABL"
                println("LLR-SPA precalculated in look-up table)")
            elseif mode == "MSUM"
                println("LLRs calculated by min-sum algorithm)")
            end
        elseif mode == "LBP"
            println("Message passing protocol: LBP")
        elseif mode == "iLBP"
            println("Message passing protocol: iLBP")
        elseif mode == "RBP"
            println("Message passing protocol: RBP")
        elseif mode == "LRBP"
            println("Message passing protocol: Local RBP")
        elseif mode == "RRBP"
            println("Message passing protocol: Randomized RBP")
        end

        println("Maximum number of iterations: $max")
        println("Simulated for SNR (dB): $snr")
        println("Stop at zero syndrome ? $stop")
        (supermode == "RBP") ? println("Decaying factor: $rbpfactor") : nothing
        (mode == "RRBP") ? println("Sample size: $SAMPLESIZE") : nothing 
        println()
    end

    if !test
        NTH = 32
        K = length(SNR)
        fer, ber = zeros(K,NTH), zeros(max,K,NTH)
        # Threads.@threads 
        for k in 1:K
            rng_noise = Xoshiro(rgn_seed_noise)
            rng_sample = (mode == "RRBP") ? Xoshiro(rng_seed_sample) : nothing
            Threads.@threads for nth in 1:NTH
                fer[k,nth], ber[:,k,nth] = 
                    performance_simulation_core(
                                        codeword,
                                        snr[k],
                                        H,
                                        mode,
                                        supermode,
                                        nreals÷NTH,
                                        max,
                                        stop,
                                        rbpfactor,
                                        rng_noise,
                                        rng_sample,
                                        test,
                                        t_test,
                                        printing)
            end
        end

        FER = zeros(K)
        BER = zeros(max,K)
        FER .= sum(fer,dims=2)/(nreals)
        BER .= sum(ber,dims=3)/(nreals*N)

        return log10.(FER), log10.(BER)
    
    else

        rng_noise = Xoshiro(rgn_seed_noise)
        rng_sample = (mode == "RRBP") ? Xoshiro(rng_seed_sample) : nothing

        Lr, Lq = performance_simulation_core(
                                    codeword,
                                    snr[1],
                                    H,
                                    mode,
                                    supermode,
                                    nreals,
                                    max,
                                    stop,
                                    rbpfactor,
                                    rng_noise,
                                    rng_sample,
                                    test,
                                    t_test,
                                    printing)
        
        return Lr, Lq
        
    end

end