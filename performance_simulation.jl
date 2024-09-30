################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER BER x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("BP.jl")
include("calc_priors.jl")
include("minsum.jl")
include("minsum_RBP.jl")

function 
    performance_simulation(
        c::Vector{Bool},
        snr::Vector{<:Real},
        H::BitMatrix,
        mode::String,
        nreals::Integer,
        max::Integer;
        t_test=nothing,
        printing=false,
        stop=false
    )

    # if nreals = 1, execute testing mode
    TEST = (nreals == 1) ? true : false
    # transform snr in standard deviations
    σ = 1 ./ sqrt.(exp10.(snr/10))

    # set random seed
    Random.seed!(SEED)

    ################################ CHECK MODE ################################
    if (mode ≠ "MKAY") && (mode ≠ "TANH") && (mode ≠ "LBP") &&
       (mode ≠ "ALTN") && (mode ≠ "TABL") && (mode ≠ "RBP") && 
       (mode ≠ "MSUM") && (mode ≠ "iLBP") && (mode ≠ "LRBP")
        throw(
            ArgumentError(
                "$mode is not a valid mode"
            )
        )
    end

    ############################### CONSTANTS ##################################

    # BPKS
    u = Float64.(2*c .- 1)
    # divisor
    divisor = nreals * N 
    # list of checks and variables nodes
    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    ############################# PREALLOCATIONS ################################

    # frame error rate
    FER = zeros(length(σ))

    # bit error rate
    ber = zeros(max)    
    BER = zeros(max,length(σ))

    # iteration in which SPA stopped
    iters = zeros(Int, length(σ), nreals)

    # MAP estimate
    d = Vector{Bool}(undef,N)

    # syndrome
    syndrome = ones(Bool,M)

    # prior llr-probabilitities
    Lf = (mode != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    # noise
    noise = Vector{Float64}(undef,N)

    # received signal
    t = Vector{Float64}(undef,N)

    # bit-error
    bit_error = Vector{Bool}(undef,N) 

    # Lq -> matrix of vn2cn messages (N x M)
    # Lr -> matrix of cn2vn messages (M x N)
    Lq, Lr = (mode != "MKAY") ? (H'*0.0,H*0.0) : (zeros(N,M,2),zeros(M,N,2))

    # Set variables that depend on the mode
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

    Ldn, visited_vns = (mode == "LBP" || mode == "iLBP" || mode == "LRBP") ?
        (zeros(N),zeros(Bool,N)) : (nothing,nothing)

    phi = (mode == "TABL") ? lookupTable() : nothing

    Residues, samples = (mode == "RBP") ? 
        (H*0.0, Vector{Int}(undef,SAMPLESIZE)) : (nothing,nothing)

    Edges, maxcoords, Factors, pfactor, num_edges  = 
        (mode == "RBP" || mode == "LRBP") ? 
        (H*0, [1,1], 1.0*H, DECAYCTE, sum(H))  : 
        (nothing,nothing,nothing,nothing,nothing)
    
    # unify the 5 flooding methods 
    _mode = (mode == "MKAY" || mode == "TANH" || mode == "ALTN" ||
             mode == "TABL" || mode == "MSUM") ? "FLOO" : mode

    ######################### FIRST RECEIVED SIGNAL ############################
    # In order to allow a test with a given received signal t_test, the first
    # received signal t is set outside the main loop.
    if TEST
        if t_test === nothing # if no test signal was provided:
            # generate a received signal
            received_signal!(t,noise,σ[1],u)
        elseif length(t_test) != N
            # if a received test signal was given but with wrong size
            throw(
                DimensionMismatch(
                    "length(t_test) should be $N, not $(length(t_test))"
                )
            )
        else
            # the received signal is the test signal
            t = t_test
        end
    else
        # generate the first received signal outside the main loop
        received_signal!(t,noise,σ[1],u)
    end

    ############################## MAIN LOOP ##################################
    if TEST 
        println(
"###################### Starting simulation(Testing mode) ######################"
    )
    else
        println(
"################## Starting simulation(# of trials: $nreals) ##################"
    )
    end
    println()
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
        print("Message passing protocol: LRBP")
    end

    println("Maximum number of iterations: $max")
    println("Simulated for SNR (dB): $snr")
    println("Stop at zero syndrome ? $stop")
    println()

    for k in eachindex(σ)
        for j in 1:nreals

            # init the llr priors
            calc_Lf!(Lf,t,σ[k]^2)
            if mode == "TABL"
                # scale for table
                Lf .*= SIZE_per_RANGE
            end            
            # initialize matrix Lr
            Lr .*= 0
            # initialize matrix Lq
            init_Lq!(Lq,Lf,vn2cn)

            if mode == "LRBP"
                # find the coordenades of the maximum residue
                minsum_lRBP_init!(maxcoords,Lq,signs,cn2vn)
            end
            if mode == "RBP"
                # initialize the matrix of residues
                minsum_RBP_init!(Residues,Lq,signs,cn2vn)
            end       
            # SPA routine
            DECODED, i = BP!(
                            _mode,
                            stop,
                            TEST,
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
                            Edges,
                            maxcoords,
                            Factors,
                            pfactor,
                            num_edges,
                            Ldn,
                            visited_vns,
                            samples
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
            received_signal!(t,noise,σ[k],u)

        end

        @inbounds @fastmath FER[k] /= NREALS

    end

    if TEST
        if mode == "RBP" || mode == "LRBP"
            return Lr, Lq, Edges
        else
            return Lr, Lq
        end
    else
        return log10.(FER), log10.(BER), iters
    end

end
