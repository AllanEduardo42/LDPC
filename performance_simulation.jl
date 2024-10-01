################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER BER x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("BP.jl")
include("calc_Lf.jl")
include("minsum.jl")
include("minsum_RBP.jl")

function
    performance_simulation(
        c::Vector{Bool},
        snr::Vector{<:Real},
        H::BitMatrix,
        mode::String,
        nreals::Integer,
        max::Integer,
        rgn_seed_noise::Integer;
        rng_seed_sample=1234,
        t_test=nothing,
        printing=false,
        stop=false
    )

########################## SET RANDOM GENERADOR SEEDS ##########################

    rng_noise = Xoshiro(rgn_seed_noise)

    rng_sample = Xoshiro(rng_seed_sample)

############################### CHECK VALID MODE ###############################
    if (mode ≠ "MKAY") && (mode ≠ "TANH") && (mode ≠ "LBP") &&
        (mode ≠ "ALTN") && (mode ≠ "TABL") && (mode ≠ "RBP") && 
        (mode ≠ "MSUM") && (mode ≠ "iLBP") && (mode ≠ "LRBP")
        throw(
            ArgumentError(
                "$mode is not a valid mode"
            )
        )
    end

################################## CONSTANTS ###################################
    
    # if nreals = 1, set test mode
    test = (nreals == 1) ? true : false
    # transform snr in standard deviations
    variances = 1 ./ (exp10.(snr/10))
    stdevs = sqrt.(variances)
    # BPKS
    u = Float64.(2*c .- 1)
    # divisor
    divisor = nreals * N 
    # list of checks and variables nodes
    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    ############################# PREALLOCATIONS ################################

    # frame error rate
    FER = zeros(length(stdevs))

    # bit error rate
    ber = zeros(max)    
    BER = zeros(max,length(stdevs))

    # iteration in which SPA stopped
    iters = Matrix{Int}(undef, length(stdevs), nreals)

    # estimate
    d = Vector{Bool}(undef,N)

    # syndrome
    syndrome = ones(Bool,M)

    # prior llr (if mode == "MKAY" just the prior probabilities)
    Lf = (mode != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    # noise
    noise = Vector{Float64}(undef,N)

    # received signal
    t = Vector{Float64}(undef,N)

    # bit-error
    bit_error = Vector{Bool}(undef,N) 

    # Lq -> matrix of vn2cn messages (N x M)
    # Lr -> matrix of cn2vn messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq, Lr = (mode != "MKAY") ? (H'*0.0,H*0.0) : (zeros(N,M,2),zeros(M,N,2))

    # Set variables Lrn and signs depending on the mode (also used for dispatch)
    if mode == "TANH" || mode == "LBP" || mode == "iLBP"
        if mode == "TANH" && !(FAST)
            Lrn = nothing
        else
            Lrn = zeros(N)
        end
        signs = nothing
    elseif mode == "ALTN" || mode == "TABL"
        Lrn = zeros(N)
        signs = zeros(Bool,N)
    elseif mode == "MSUM" || mode == "RBP" || mode == "LRBP"
        Lrn = nothing
        signs = zeros(Bool,N)
    else
        Lrn = nothing
        signs = nothing
    end

   # Set other variables that depend on the mode

    visited_vns = (mode == "LBP" || mode == "iLBP") ? zeros(Bool,N) : nothing

    Ldn = (mode == "LBP" || mode == "iLBP" || mode == "RBP" || mode == "LRBP") ?
          zeros(N) : nothing

    phi = (mode == "TABL") ? lookupTable() : nothing

    Residues = (mode == "RBP") ? H*0.0 : nothing

    samples = (mode == "RBP" && SAMPLESIZE != 0) ?
                Vector{Int}(undef,SAMPLESIZE) : nothing

    if mode == "RBP"
        rbpfactor = DECAYRBP
    elseif mode == "LRBP"
        rbpfactor = DECAYLRBP
    else
        rbpfactor = nothing
    end

    maxcoords, Factors, num_edges  = 
        (mode == "RBP" || mode == "LRBP") ? 
        ([1,1], 1.0*H, sum(H)) : 
        (nothing,nothing,nothing)

    # unify the 5 flooding methods 
    _mode = (mode == "MKAY" || mode == "TANH" || mode == "ALTN" ||
                mode == "TABL" || mode == "MSUM") ? "FLOO" : mode

    ######################### FIRST RECEIVED SIGNAL ############################
    # In order to allow to test with a given received signal t_test, the first
    # received signal t must be set outside the main loop.
    if test
        if t_test === nothing # if no test signal was provided:
            # generate a new received signal
            received_signal!(t,noise,stdevs[1],u,rng_noise)
        elseif length(t_test) != N
            # if a received test signal was given, but with wrong length
            throw(
                DimensionMismatch(
                    "length(t_test) should be $N, not $(length(t_test))"
                )
            )
        else
            # the received signal is the given test signal
            t = t_test
        end
    else
        # generate the first received signal outside the main loop
        received_signal!(t,noise,stdevs[1],u,rng_noise)
    end

########################## PRINT SIMULATION DETAILS ############################
    if test && printing
        println()
        println("""###################### Starting simulation (Testing mode) ###
        ###################""")
        println()
    elseif !test
        println()
        println("""############################# Starting simulation ###########
        ###################""")
        println()
        println("Number of trials: $nreals")
    end
    if !test || printing
        if _mode == "FLOO"
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
            println("Message passing protocol: LRBP")
        end

        println("Maximum number of iterations: $max")
        println("Simulated for SNR (dB): $snr")
        println("Stop at zero syndrome ? $stop")
        println()
    end

################################## MAIN LOOP ###################################
    
    for k in eachindex(stdevs)
        for j in 1:nreals

            # init the llr priors
            calc_Lf!(Lf,t,variances[k])
            if mode == "TABL"
                # scale for table
                Lf .*= SIZE_per_RANGE
            end            
            # initialize matrix Lr
            Lr .*= 0
            # initialize matrix Lq
            init_Lq!(Lq,Lf,vn2cn)

            if mode == "RBP" || mode == "LRBP"
                minsum_RBP_init!(Residues,maxcoords,Lq,signs,cn2vn)
            end      
            # SPA routine
            DECODED, i = BP!(
                            _mode,
                            stop,
                            test,
                            max,
                            syndrome,
                            d,
                            c,
                            bit_error,
                            ber,
                            Lf,
                            Lq,
                            Lr,
                            cn2vn,
                            vn2cn,
                            Lrn,
                            signs,
                            phi,
                            printing,
                            Residues,
                            maxcoords,
                            Factors,
                            rbpfactor,
                            num_edges,
                            Ldn,
                            visited_vns,
                            samples,
                            rng_sample
                        )                

            # bit error rate
            @fastmath ber ./= divisor
            @inbounds  @fastmath @. BER[:,k] += ber
            # iteration in which SPA stopped (iszero(syndrome) = true)
            @inbounds iters[k,j] = i
            if !(DECODED)
                # frame error rate
                @inbounds FER[k] += 1
            end

            # received signal for the next realization (j+1)
            received_signal!(t,noise,stdevs[k],u,rng_noise)

        end

        @inbounds @fastmath FER[k] /= NREALS

    end

    if test
        return Lr, Lq
    else
        return log10.(FER), log10.(BER), iters
    end

end
