################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER bit_error SNR)

function 
    performance_estimation(
        c::Vector{Bool},
        σ::Vector{<:AbstractFloat},
        H::BitMatrix,
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}, 
        phi::Vector{<:AbstractFloat},
        mode::String;
        nreals = NREALS
    )

    ############################### constants ##################################

    M = length(checks2nodes)
    N = length(nodes2checks)
    divisor = NREALS * N
    # BPKS
    u = Float64.(2*c .- 1)

    ############################# preallocation ################################

    # frame error rate
    FER = zeros(length(σ))

    # bit error rate
    ber = zeros(MAX)
    BER = zeros(MAX,length(σ))

    # iteration in which SPA stopped
    iters = zeros(Int, length(σ), NREALS)

    # prior Δ-llr
    ΔLf = Vector{Float64}(undef,N)

    # Vertical and horizontal update matrices
    Lq = H*0.0
    Lr = H*0.0
    Lr_return = H*0.0

    # received signal
    t = Vector{Float64}(undef,N)

    # MAP estimate
    d = Vector{Bool}(undef,N)

    # syndrome
    syndrome = Vector{Bool}(undef,M)

    # noise
    noise = Vector{Float64}(undef,N)

    # auxiliary variables
    Lrn = zeros(N)
    sn = ones(Int8,N)
    bit_error = Vector{Bool}(undef,N)    

    ################################## MAIN ####################################

    for k in eachindex(σ)
        
        # set random seed
        Random.seed!(SEED)

        σ² = σ[k]^2
    
        for j in 1:nreals

            randn!(noise)

            noise .*= σ[k]

            t .= u .+ noise

            # prior Δ-llr

            for i in eachindex(t)
                if mode == "TAB"
                    @inbounds ΔLf[i] = -2*SIZE_per_RANGE*t[i]/σ²
                else
                    @inbounds ΔLf[i] = -2t[i]/σ²
                end
            end
            
            # initialize matrix Lq
            llr_init_q!(Lq,ΔLf,nodes2checks)
            
            # SPA
            i = 0
            DECODED = false
            if mode == "TNH"
                # tanh SPA
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
                        ΔLf,
                        syndrome,
                        nothing,
                        Lrn,
                        nothing
                    )
                ;
            elseif mode == "ALT"
                # alternative SPA
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
                        ΔLf,
                        syndrome,
                        sn,
                        Lrn,
                        nothing
                    )
                ;
            elseif mode == "TAB"
                # lookup-table SPA
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
                        ΔLf,
                        syndrome,
                        sn,
                        Lrn,
                        phi
                    )
                ;
            elseif mode == "MIN"
                # approximate SPA
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
                        ΔLf,
                        syndrome,
                        sn,
                        nothing,
                        nothing
                    )
                ;
            elseif mode == "LBP"
                # LBP SPA
                Lr = zeros(M,N)
                DECODED, i = 
                    SPA_LBP!(
                        d,
                        ber,
                        c,
                        bit_error,
                        Lr,
                        Lq,
                        checks2nodes,
                        nodes2checks,
                        ΔLf,
                        syndrome,
                        Lrn
                    )
                ;
            elseif mode == "RBP"
                # RBP SPA
                Lr = 0.0*H
                DECODED, i = 
                    SPA_RBP!(
                        d,
                        ber,
                        c,
                        bit_error,
                        Lr,
                        Lq,
                        checks2nodes,
                        nodes2checks,
                        ΔLf,
                        syndrome,
                        Lrn
                    )
                ;
            else
                throw(
                    ArgumentError(
                        "$mode is not a valid mode"
                    )
                )
            end
            # bit error rate
            @fastmath ber ./= divisor
            @inbounds  @fastmath @. BER[:,k] += ber
            # iteration in which SPA stopped
            @inbounds iters[k,j] = i
            if ~DECODED
                # frame error rate
                @inbounds FER[k] += 1
            end
        end

        @inbounds @fastmath FER[k] /= NREALS

    end

    if nreals == 1
        return Lr, Lq
    else
        return log10.(FER), log10.(BER), iters
    end
end