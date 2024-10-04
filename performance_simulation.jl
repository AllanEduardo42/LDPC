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

        K = length(SNR)
        FER, BER, Iters = zeros(K), zeros(max,K), zeros(K,nreals)
        Threads.@threads for i in eachindex(SNR)
            FER[i], BER[:,i], Iters[i,:] = 
                performance_simulation_core(
                                    codeword,
                                    snr[i],
                                    H,
                                    mode,
                                    supermode,
                                    nreals,
                                    max,
                                    stop,
                                    rbpfactor,
                                    rgn_seed_noise,
                                    rng_seed_sample,
                                    test,
                                    t_test,
                                    printing)
        end

        return FER, BER, Iters
    
    else

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
                                    rgn_seed_noise,
                                    rng_seed_sample,
                                    test,
                                    t_test,
                                    printing)
        
        return Lr, Lq
        
    end

end