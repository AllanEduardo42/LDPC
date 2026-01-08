################################################################################
# Allan Eduardo Feitosa
# 26 mai 2024
# First layer of the routine to estimate the LPDC performance (Fer Ber x SNR)

include("simcore.jl")
include("print_simulation_details.jl")

function 
    prepare_simulation(
        ebn0::Vector{Float64},
        algorithm::String,
        trials::Vector{Int},
        maxiter::Integer,
        bptype::String,
        decay::Float64
    )

    if bptype == "MKAY" && algorithm ≠ "Flooding"
        throw(error(lazy"The McKay method is only supported for Flooding"))
    end

########################### PRINT SIMULATION DETAILS ###########################
    print_simulation_details(TEST,trials,algorithm,bptype,maxiter,ebn0,decay)

################################ MULTITHREADING ################################   
    Lr = Matrix{Float64}(undef,MM,NN)
    Lq = Matrix{Float64}(undef,MM,NN)
    if TEST
        K = 1
        nthreads = 1
    else
        K = length(ebn0) 
        nthreads = NTHREADS
    end
    sum_greediness = Array{Int,4}(undef,maxiter,sum(HH)+1,K,nthreads)
    sum_decoded = Array{Int,3}(undef,maxiter,K,nthreads)
    sum_ber = Array{Int,3}(undef,maxiter,K,nthreads)
    for k in 1:K
        stats = @timed Threads.@threads for i in 1:nthreads
        # stats = @timed for i in 1:nthreads
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
                                algorithm,
                                bptype,
                                trials[k]÷NTHREADS,
                                maxiter,
                                STOP,
                                decay,
                                LISTSIZES,
                                RGN_SEEDS[i],
                                TEST,
                                PRIN
                            )
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

        prob_greediness = zeros(maxiter,sum(HH)+1,K)
        prob_greediness .= sum(sum_greediness,dims=4)
        for k = 1:K
            prob_greediness[:,:,k] ./= trials[k]
        end

        println()

        lowerfer = 1/maximum(trials)
        lowerber = lowerfer/NN
        replace!(x-> x < lowerfer ? lowerfer : x, Fer)
        replace!(x-> x < lowerber ? lowerber : x, Ber)

        return Fer, Ber, prob_greediness
    end
end