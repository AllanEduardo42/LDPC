################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER BER x SNR)

include("lookupTable.jl")
include("SPA.jl")
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
    if (mode ≠ "TNH") && (mode ≠ "ALT") && (mode ≠ "TAB") && (mode ≠ "MIN") &&
        (mode ≠ "LBP") && (mode ≠ "RBP") && (mode ≠ "RBP_R")
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
    syndrome = Vector{Bool}(undef,M)

    # prior llr-probabilitities
    Lf = Vector{Float64}(undef,N)

    # noise
    noise = Vector{Float64}(undef,N)

    # received signal
    t = Vector{Float64}(undef,N)

    # bit-error
    bit_error = Vector{Bool}(undef,N) 

    # Vertical and horizontal update matrices
    Lq = H*0.0
    Lr = H*0.0

    # variables used only if testing
    f = TEST ? zeros(N,2) : nothing
    q,r = TEST ? (zeros(M,N,2),zeros(M,N,2)) : (nothing,nothing)

    # Set variables that depend on the mode
    if mode == "TNH" || mode == "LBP"
        Lrn = zeros(N)
        sn = nothing
    elseif mode == "ALT" || mode == "TAB"
        Lrn = zeros(N)
        sn = zeros(Bool,N)
    elseif mode == "MIN" || mode == "RBP" || mode == "RBP_R"
        Lrn = nothing
        sn = zeros(Bool,N)
    end

    phi = (mode == "TAB") ? lookupTable() : nothing
    R = (mode == "RBP_R") ? H*0.0 : nothing
    Edges = (mode == "RBP" || mode == "RBP_R") ? H*0 : nothing
    max_coords = (mode == "RBP" || mode == "RBP_R") ? [1,1] : nothing
    penalty = (mode == "RBP" || mode == "RBP_R") ? 1.0*H : nothing
    penalty_factor = (mode == "RBP" || mode == "RBP_R") ? PENALTY : nothing
    num_edges = (mode == "RBP" || mode == "RBP_R") ? sum(H) : nothing

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
        # init the priors for test
        calc_f!(f,t,σ[1]^2)
        # init matrix q
        init_q!(q,f,nodes2checks)
        # init matrix r
        r = zeros(M,N,2)
    else
        # generate the first received signal outside the main loop
        received!(t,noise,σ[1],u)
    end

    ############################## MAIN LOOP ##################################
    for k in eachindex(σ)
        for j in 1:nreals

            # init the llr priors
            calc_Lf!(Lf,t,σ[k]^2)
            if mode == "TAB"
                # scale for table
                Lf .*= SIZE_per_RANGE
            end            
            # initialize matrix Lr
            Lr .*= 0
            # initialize matrix Lq
            llr_init_q!(Lq,Lf,nodes2checks)

            if mode == "RBP"
                # find max_coords for the first update
                min_sum_RBP_init!(
                    max_coords,
                    Lq,
                    sn,
                    checks2nodes
                )
            end
            if mode == "RBP_R"
                # initialize the matrix of residues R
                min_sum_RBP_R_init!(
                    R,
                    Lq,
                    sn,
                    checks2nodes
                )
            end       
            # SPA routine
            DECODED, i = 
            SPA!(
                mode,
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
                f,
                q,
                r,
                printing,
                R,
                Edges,
                max_coords,
                penalty,
                penalty_factor,
                num_edges
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
        if mode == "RBP" || mode == "RBP_R"
            return r, Lr, q, Lq, Edges
        else
            return r, Lr, q, Lq
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
