################################################################################
# Allan Eduardo Feitosa
# 26 mai 2024
# First layer of the routine to estimate the LPDC performance (Fer Ber x EbN0)
# Separates simulations in terms of EbN0 and number of threads (paralelization)
# Returns FER and BER matrices (iteration x EbN)

function prepare_simulation(
    count::Int,
    ebn0::Vector{Float64},
    algorithm::String,
    max_frame_errors::Int,
    maxiter::Integer,
    decay::Float64
)

############################ PRINT ALGORITHM DETAILS ###########################
    print_algorithm_details(count,algorithm,BPTYPE,decay)

################################ MULTITHREADING ################################   
    if TEST
        Lr = Matrix{Float64}(undef,MM,NN)
        Lq = Matrix{Float64}(undef,MM,NN)
        num_ebn0 = 1
        nthreads = 1
    else
        Lr = nothing
        Lq = nothing
        num_ebn0 = length(ebn0) 
        nthreads = NTHREADS
    end
    max_frame_errors ÷= nthreads
    Frame_errors = Array{Int,3}(undef,maxiter,num_ebn0,nthreads)
    Bit_errors = Array{Int,3}(undef,maxiter,num_ebn0,nthreads)
    Trials = Matrix{Int}(undef,num_ebn0,nthreads)
    for k in 1:num_ebn0

        # transform EbN0 in standard deviations
        variance = exp10.(-ebn0[k]/10) / (2*RR)
        stdev = sqrt.(variance)

        # thread "for" loop
        stats = @timed Threads.@threads for i in 1:nthreads
        # stats = @timed for i in 1:nthreads
            Lr, Lq, frame_errors, bit_errors, trials = simcore(
                                                    KK,
                                                    CODE_LENGTH,
                                                    G_CRC,
                                                    stdev,
                                                    HH,
                                                    PP,
                                                    NC,
                                                    NV,
                                                    E_H,
                                                    PROTOCOL,
                                                    LIFTSIZE,
                                                    algorithm,
                                                    BPTYPE,
                                                    max_frame_errors,
                                                    maxiter,
                                                    RAYL,
                                                    C_DR_ITER,
                                                    decay,
                                                    LISTSIZES,
                                                    CI_GAMMA,
                                                    RGN_SEEDS[i],
                                                    TEST,
                                                    TEST && PRIN
                                                )
            if !TEST
                Frame_errors[:,k,i] = frame_errors
                Bit_errors[:,k,i] = bit_errors
                Trials[k,i] = trials
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
        Total_trials = sum(Trials,dims=2)
        Fer = zeros(maxiter,num_ebn0)
        Fer .= sum(Frame_errors,dims=3)
        for k = 1:num_ebn0
            Fer[:,k] ./= Total_trials[k]
        end
        Ber = zeros(maxiter,num_ebn0)
        Ber .= sum(Bit_errors,dims=3)
        for k = 1:num_ebn0
            Ber[:,k] ./= (NN*Total_trials[k])
        end

        println()

        lowerfer = 1/maximum(maximum(Total_trials))
        lowerber = lowerfer/NN
        replace!(x-> x < lowerfer ? lowerfer : x, Fer)
        replace!(x-> x < lowerber ? lowerber : x, Ber)

        return Fer, Ber
    end
end