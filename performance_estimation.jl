################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER BER x SNR)

include("lookupTable.jl")
include("BP.jl")
include("calc_priors.jl")
include("min_sum.jl")
include("min_sum_RBP.jl")

function 
    performance_estimation(
        c::Vector{Bool},
        σ::Vector{<:AbstractFloat},
        H::BitMatrix,
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
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
    if (mode ≠ "MKAY") && (mode ≠ "FTNH") && (mode ≠ "FALT") &&
       (mode ≠ "FTAB") && (mode ≠ "FMSM") && (mode ≠ "oLBP") &&
       (mode ≠ "iLBP") && (mode ≠ "oRBP") && (mode ≠ "LRBP")
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

    # Vertical and horizontal update matrices
    Lq, Lr = (mode != "MKAY") ? (H*0.0,H*0.0) : (zeros(M,N,2),zeros(M,N,2))

    # Set variables that depend on the mode
    if mode == "FTNH" || mode == "oLBP" || mode == "iLBP"
        Lrn = zeros(N)
        sn = nothing
    elseif mode == "FALT" || mode == "FTAB"
        Lrn = zeros(N)
        sn = zeros(Bool,N)
    elseif mode == "FMSM" || mode == "oRBP" || mode == "LRBP"
        Lrn = nothing
        sn = zeros(Bool,N)
    else
        Lrn = nothing
        sn = nothing
    end
 
    Ldn, visited_nodes = (mode == "oLBP" || mode == "iLBP") ?
        (zeros(N),zeros(Bool,N)) : (nothing,nothing)

    phi = (mode == "FTAB") ? lookupTable() : nothing

    Residues = (mode == "oRBP") ? H*0.0 : nothing

    Edges, maxcoords, Factors, pfactor, num_edges  = 
        (mode == "oRBP" || mode == "LRBP") ? 
        (H*0, [1,1], 1.0*H, PENALTY, sum(H))  : 
        (nothing,nothing,nothing,nothing,nothing)
    
    # unity the 5 flooding methods 
    _mode = (mode == "MKAY" || mode == "FTNH" || mode == "FALT" ||
             mode == "FTAB" || mode == "FMSM") ? "FLOO" : mode

    ######################### FIRST RECEIVED SIGNAL ############################
    # In order to allow a test with a given received signal t_test, the first
    # received signal t is set outside the main loop.
    if TEST
        if t_test === nothing # if no test signal was provided:
            # generate a received signal
            received!(t,noise,σ[1],u)
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
        received!(t,noise,σ[1],u)
    end

    ############################## MAIN LOOP ##################################
    for k in eachindex(σ)
        for j in 1:nreals

            # init the llr priors
            calc_Lf!(Lf,t,σ[k]^2)
            if mode == "FTAB"
                # scale for table
                Lf .*= SIZE_per_RANGE
            end            
            # initialize matrix Lr
            Lr .*= 0
            # initialize matrix Lq
            init_Lq!(Lq,Lf,nodes2checks)

            if mode == "LRBP"
                # find the coordenades of the maximum residue
                min_sum_lRBP_init!(
                    maxcoords,
                    Lq,
                    sn,
                    checks2nodes
                )
            end
            if mode == "oRBP"
                # initialize the matrix of residues
                min_sum_RBP_init!(
                    Residues,
                    Lq,
                    sn,
                    checks2nodes
                )
            end       
            # SPA routine
            DECODED, i = 
            BP!(
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
                checks2nodes,
                nodes2checks,
                Lrn,
                sn,
                phi,
                printing,
                Residues,
                Edges,
                maxcoords,
                Factors,
                pfactor,
                num_edges,
                Ldn,
                visited_nodes
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
            received!(t,noise,σ[k],u)

        end

        @inbounds @fastmath FER[k] /= NREALS

    end

    if TEST
        if mode == "oRBP" || mode == "LRBP"
            return Lr, Lq, Edges
        else
            return Lr, Lq
        end
    else
        return log10.(FER), log10.(BER), iters
    end

end

function
    received!(
        t::Vector{<:AbstractFloat},
        noise::Vector{<:AbstractFloat},
        σ::AbstractFloat,
        u::Vector{<:AbstractFloat})

    randn!(noise)
    @fastmath noise .*= σ
    @fastmath t .= u .+ noise

end
