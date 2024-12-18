################################################################################
# Allan Eduardo Feitosa
# 3 out 2024
# First layer of the routine to estimate the LPDC performance (FER BER x SNR)

include("simcore.jl")

function 
    performance_sim(
        snr::Vector{<:Real},
        mode::String,
        trials::Vector{<:Integer},
        maxiter::Integer,
        bptype::String,
        decayfactor::AbstractFloat;
        compile = false,
        test = false,
        printtest = false    
    )

########################### PRINT SIMULATION DETAILS ###########################
    if !compile || test
        println()
        print("############################### LDPC parameters #######################")
        println("#########")
        println()
        println("Parity Check Matrix: $M x $N")
        println()
        display(sparse(H))
        println()
        println("Graph girth = ", girth)
        println()
    end
    if test
        print("###################### Starting simulation (Testing mode) #####")
        println("#################")
        println()
    elseif !compile
        print("############################# Starting simulation #############")
        println("#################")
        println()
    end
    if !compile || test
        println("Number of trials: $trials")
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
        println("Number of threads (multithreading): $NTHREADS")
        println("Simulated for SNR (dB): $snr")
        println("Stop at zero syndrome ? $STOP")
        (mode == "RBP") || (mode == "Local-RBP") || (mode == "List-RBP") ?
        println("RBP decaying factor: $decayfactor") : nothing
        (mode == "List-RBP") ? println("List 1 size: $LISTSIZE\nList 2 size: $LISTSIZE2") : nothing 
        println()
    end

################################ MULTITHREADING ################################

    if compile
        Lr, Lq = performance_simcore(
            A,
            snr[1],
            H,
            E_H,
            zf,
            Zc,
            mode,
            bptype,
            trials[1],
            maxiter,
            STOP,
            decayfactor,
            LISTSIZE,
            LISTSIZE2,
            Rgn_noise_seeds[1],
            Rgn_samples_seeds[1],
            Rgn_message_seeds[1];
            test=true,
            printtest=printtest)

        return Lr, Lq 
    else
        K = length(SNR)
        if MTHR
            decoded, ber = zeros(maxiter,K,NTHREADS), zeros(maxiter,K,NTHREADS)
        else
            decoded, ber = zeros(maxiter,K), zeros(maxiter,K)
        end
        for k in 1:K
            if MTHR
                @time Threads.@threads for i in 1:NTHREADS
                    decoded[:,k,i], ber[:,k,i] = performance_simcore(
                                            A,
                                            snr[k],
                                            H,
                                            E_H,
                                            zf,
                                            Zc,
                                            mode,
                                            bptype,
                                            trials[k]÷NTHREADS,
                                            maxiter,
                                            STOP,
                                            decayfactor,
                                            LISTSIZE,
                                            LISTSIZE2,
                                            Rgn_noise_seeds[i],
                                            Rgn_samples_seeds[i],
                                            Rgn_message_seeds[i])
                end
            else
                @time decoded[:,k], ber[:,k] = performance_simcore(
                                            A,
                                            snr[k],
                                            H,
                                            E_H,
                                            zf,
                                            Zc,
                                            mode,
                                            bptype,
                                            trials[k],
                                            maxiter,
                                            STOP,
                                            decayfactor,
                                            LISTSIZE,
                                            LISTSIZE2,
                                            Rgn_noise_seeds[1],
                                            Rgn_samples_seeds[1],
                                            Rgn_message_seeds[1])
            end
        end

        if MTHR
            FER = sum(decoded,dims=3)
            for k = 1:K
                FER[:,k] ./= trials[k]
            end
            @. FER = 1 - FER
            BER = sum(ber,dims=3)
            for k = 1:K
                BER[:,k] ./= (N*trials[k])
            end
        else
            FER = similar(decoded)
            BER = similar(ber)
            for k = 1:K
                FER[:,k] = decoded[:,k]/trials[k]
            end
            FER = 1 .- FER
            for k = 1:K
                BER[:,k] = ber[:,k]/(N*trials[k])
            end
        end

        return log10.(FER), log10.(BER)
       
    end
end