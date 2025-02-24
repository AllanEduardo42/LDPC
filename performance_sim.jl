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
        decay::AbstractFloat,
        stop::Bool,
        mthr::Bool;
        test = false,
        printtest = false    
    )

########################### PRINT SIMULATION DETAILS ###########################
    if test
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
        println(file,str)
    end
    str = "Number of trials: $trials"
    println(str)
    if SAVE
        println(file,str)
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
        println(file,str)
    end
    str=
    """Maximum number of iterations: $maxiter
    Number of threads (multithreading): $NTHREADS
    Simulated for SNR (dB): $snr
    Stop at zero syndrome ? $stop"""
    println(str)
    if SAVE
        println(file,str)
    end
    if mode == "RBP" || mode == "Local-RBP" || mode == "List-RBP" ||
       mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        str = "RBP decaying factor: $decay"
        println(str)
        if SAVE
            println(file,str)
        end
    end
    if mode == "List-RBP" || mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        str = "List 1 size: $(Listsizes[1])\nList 2 size: $(Listsizes[2])"
        println(str)
        if SAVE
            println(file,str)
        end
    end
    println()
    if SAVE
        println(file,"")
    end

    ################################ MULTITHREADING ################################

    if test
        Lr, Lq, max_residues = simcore(
            A,
            snr[1],
            H,
            E_H,
            LDPC,
            Zf,
            nr_ldpc_data,
            mode,
            bptype,
            trials[1],
            maxiter,
            stop,
            decay,
            Listsizes,
            Rgn_noise_seeds[1],
            Rgn_samples_seeds[1],
            Rgn_message_seeds[1];
            test=true,
            printtest=printtest)
        
        println()
        # aux = filter(isfinite,max_residues)
        # maxi = maximum(aux)
        # replace!(x -> isfinite(x) ? x : 2maxi, max_residues)

        return Lr, Lq, max_residues

    else
        K = length(snr)
        if mthr
            decoded, ber = zeros(maxiter,K,NTHREADS), zeros(maxiter,K,NTHREADS)
        else
            decoded, ber = zeros(maxiter,K), zeros(maxiter,K)
        end
        for k in 1:K
            if mthr
                stats = @timed Threads.@threads for i in 1:NTHREADS
                    decoded[:,k,i], ber[:,k,i] = simcore(
                                                    A,
                                                    snr[k],
                                                    H,
                                                    E_H,
                                                    LDPC,
                                                    Zf,
                                                    nr_ldpc_data,
                                                    mode,
                                                    bptype,
                                                    trials[k]÷NTHREADS,
                                                    maxiter,
                                                    stop,
                                                    decay,
                                                    Listsizes,
                                                    Rgn_noise_seeds[i],
                                                    Rgn_samples_seeds[i],
                                                    Rgn_message_seeds[i])
                end
            else
                stats = @timed decoded[:,k], ber[:,k] = simcore(
                                                    A,
                                                    snr[k],
                                                    H,
                                                    E_H,
                                                    LDPC,
                                                    Zf,
                                                    nr_ldpc_data,
                                                    mode,
                                                    bptype,
                                                    trials[k],
                                                    maxiter,
                                                    stop,
                                                    decay,
                                                    Listsizes,
                                                    Rgn_noise_seeds[1],
                                                    Rgn_samples_seeds[1],
                                                    Rgn_message_seeds[1])
            end
            str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
            println(str)
            if SAVE
                println(file,str)
            end
        end

        if mthr
            FER = sum(decoded,dims=3)[:,:,1]
            for k = 1:K
                FER[:,k] ./= trials[k]
            end
            @. FER = 1 - FER
            BER = sum(ber,dims=3)[:,:,1]
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

        println()

        lower = 1/maximum(trials)
        replace!(x-> x < lower ? lower : x, FER)
        replace!(x-> x < lower ? lower : x, BER)

        return log10.(FER), log10.(BER)
        
    end

end