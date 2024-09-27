################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER BER x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("BP.jl")
include("calc_priors.jl")
include("min_sum.jl")
include("min_sum_RBP.jl")

function 
    performance_simulation(
        c::Vector{Bool},
        σ::Vector{<:AbstractFloat},
        H::BitMatrix,
        mode::String,
        nreals::Integer,
        max::Integer;
        t_test=nothing,
        printing=false,
    )

    TEST = (nreals == 1) ? true : false

    # set random seed
    Random.seed!(SEED)

    ################################ CHECK MODE ####################################
    if (mode ≠ "MKAY") && (mode ≠ "TANH") && (mode ≠ "LBP") &&
       (mode ≠ "ALTN") && (mode ≠ "TABL") && (mode ≠ "RBP") && 
       (mode ≠ "MSUM") && (mode ≠ "iLBP") && (mode ≠ "LRBP")
        throw(
            ArgumentError(
                "$mode is not a valid mode"
            )
        )
    end

    ############################### constants ##################################

    # BPKS
    u = Float64.(2*c .- 1)
    # divisor
    divisor = nreals * N

    ############################# AUXILIARY CONSTANTS ##############################

    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    ############################# preallocation ################################

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
        Lrn = zeros(N)
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
 
    Ldn, visited_vns = (mode == "LBP" || mode == "iLBP") ?
        (zeros(N),zeros(Bool,N)) : (nothing,nothing)

    phi = (mode == "TABL") ? lookupTable() : nothing

    Residues = (mode == "RBP") ? H*0.0 : nothing

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
                min_sum_lRBP_init!(maxcoords,Lq,signs,cn2vn)
            end
            if mode == "RBP"
                # initialize the matrix of residues
                min_sum_RBP_init!(Residues,Lq,signs,cn2vn)
            end       
            # SPA routine
            DECODED, i = BP!(
                            _mode,
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
                            visited_vns
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
