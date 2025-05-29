################################################################################
# Allan Eduardo Feitosa
# 26 mai 2024
# First layer of the routine to estimate the LPDC performance (Fer Ber x SNR)

include("simcore.jl")
include("print_simulation_details.jl")

function 
    prepare_simulation(
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
    print_simulation_details(TEST,trials,mode,bptype,maxiter,ebn0,decay)

################################ MULTITHREADING ################################   
    if TEST
        Lr = zeros(MM,NN)
        Lq = zeros(MM,NN)
        K = 1
        nthreads = 1
    else
        K = length(ebn0) 
        nthreads = NTHREADS
        sum_decoded = zeros(Int,maxiter,K,nthreads)
        sum_ber = zeros(Int,maxiter,K,nthreads)
    end
    for k in 1:K
        stats = @timed Threads.@threads for i in 1:nthreads
            x, y, z, w = simcore(
                                AA,
                                KK,
                                RR,
                                GG,
                                G_CRC,
                                ebn0[k],
                                HH,
                                H1,
                                LL,
                                UU,
                                NC,
                                NV,
                                E_H,
                                PROTOCOL,
                                LIFTSIZE,
                                mode,
                                bptype,
                                trials[k]÷NTHREADS,
                                maxiter,
                                STOP,
                                decay,
                                LISTSIZES,
                                RELATIVE,
                                RGN_NOISE_SEEDS[i],
                                RGN_MESSAGE_SEEDS[i],
                                TEST,
                                PRIN,
                                MSGTEST,
                                NOISETEST)
            if TEST
                Lr .= x
                Lq .= y
            else
                sum_decoded[:,k,i] = z
                sum_ber[:,k,i] = w
            end
        end
        str = """Elapsed $(round(stats.time;digits=1)) seconds ($(round(stats.gctime/stats.time*100;digits=2))% gc time, $(round(stats.compile_time/stats.time*100,digits=2))% compilation time)"""
        println(str)
        if SAVE
            println(FILE,str)
        end
    end
    if TEST
        return Lr, Lq
    else
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
    end
end