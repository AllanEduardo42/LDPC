################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER BER x SNR)

include("lookupTable.jl")
include("SPA.jl")
include("calc_priors.jl")

function 
    performance_estimation(
        c::Vector{Bool},
        σ::Vector{<:AbstractFloat},
        H::BitMatrix,
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        mode::String,
        nreals::Integer;
        test=false,
        t_test=nothing,
        printing=nothing,
    )

    # set random seed
    Random.seed!(SEED)

    ############################## CHECK MODE ##################################
    if (mode ≠ "TNH") && (mode ≠ "ALT") && (mode ≠ "TAB") && (mode ≠ "MIN") &&
       (mode ≠ "LBP") && (mode ≠ "RBP")
        throw(
            ArgumentError(
                "$mode is not a valid mode"
            )
        )
    end

    ############################### constants ##################################

    M = length(checks2nodes)
    N = length(nodes2checks)
    # BPKS
    u = Float64.(2*c .- 1)

    ############################# preallocation ################################

    # constant
    divisor = nreals * N

    # frame error rate
    FER = zeros(length(σ))

    # bit error rate
    ber = zeros(MAX)    
    BER = zeros(MAX,length(σ))

    # iteration in which SPA stopped
    iters = zeros(Int, length(σ), nreals)

    # MAP estimate
    d = Vector{Bool}(undef,N)
    d_test = (test ? Vector{Bool}(undef,N) : nothing)

    # syndrome
    syndrome = Vector{Bool}(undef,M)
    syndrome_test = (test ? Vector{Bool}(undef,M) : nothing)

    # prior llr-probabilitities
    Lf = Vector{Float64}(undef,N)
    f = (test ? Matrix{Float64}(undef,N,2) : nothing)

    # noise
    noise = Vector{Float64}(undef,N)

    # received signal
    t = Vector{Float64}(undef,N)

    # bit-error
    bit_error = Vector{Bool}(undef,N) 

    # Vertical and horizontal update matrices
    Lq = H*0.0
    Lr = H*0.0
    r, q = (test ? (zeros(M,N,2), zeros(M,N,2)) : (nothing, nothing))
   
    # Set variables that depend on the mode chosen
    if mode == "TNH"
        flooding = true
        Lrn = zeros(N)
        sn = nothing
        phi = nothing
    elseif mode == "ALT"
        flooding = true
        Lrn = zeros(N)
        sn = ones(Int8,N)
        phi = nothing
    elseif mode == "TAB"
        flooding = true
        Lrn = zeros(N)
        sn = ones(Int8,N)
        phi = lookupTable()
    elseif mode == "MIN"
        flooding = true
        Lrn = nothing
        sn = ones(Int8,N)
        phi = nothing
    elseif mode == "LBP" || mode == "RBP"
        flooding = false
        Lrn = zeros(N)
        sn = nothing
        phi = nothing
    end
    
    # first realization
    if test
        if t_test === nothing
            received!(t,noise,σ[1],u)
        elseif length(t_test) != N
            throw(
                DimensionMismatch(
                    "length(t_test) should be $N, not $(length(t_test))"
                )
            )
        else
            t = t_test
        end
        calc_f!(f,t,σ[1]^2)
        init_q!(q,f,nodes2checks)
        nreals = 1
    else
        received!(t,noise,σ[1],u)
    end

    ################################## MAIN ####################################
    for k in eachindex(σ)    
        for j in 1:nreals
            # prior f-llr
            calc_Lf!(
                Lf,
                t,
                σ[k]^2
            )
            if mode == "TAB"
                Lf .*= SIZE_per_RANGE
            end
            
            # initialize matrix Lr
            Lr .*= 0
            # initialize matrix Lq
            llr_init_q!(Lq,Lf,nodes2checks)                
            # SPA
            DECODED, i = 
                SPA!(
                    d,
                    ber,
                    c,
                    bit_error,
                    Lr,
                    Lq,
                    checks2nodes,
                    nodes2checks,
                    Lf,
                    syndrome,
                    Lrn,
                    sn,
                    phi,
                    mode,
                    flooding,
                    test,
                    d_test,
                    syndrome_test,
                    f,
                    r,
                    q,
                    printing
                )                

            # bit error rate
            @fastmath ber ./= divisor
            @inbounds  @fastmath @. BER[:,k] += ber
            # iteration in which SPA stopped
            @inbounds iters[k,j] = i
            if ~DECODED
                # frame error rate
                @inbounds FER[k] += 1
            end
            # next realization
            received!(t,noise,σ[k],u)
        end

        @inbounds @fastmath FER[k] /= NREALS

    end

    if test
        return r, Lr, q, Lq
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
