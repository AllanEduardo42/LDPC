################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER bit_error SNR)

function 
    performance_estimation(
        c::Vector{Bool},
        σ::Vector{Float64},
        H::BitMatrix,
        indices_row::Vector{Vector{Int64}},
        indices_col::Vector{Vector{Int64}}, 
        phi::Vector{Float64},
        mode::String;
        nreals = NREALS
    )

    ############################### constants ##################################

    M = length(indices_row)
    N = length(indices_col)
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
    Lq = H'*0.0
    Lr = H*0.0

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
    sn = ones(Bool,N)
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
            llr_init_q!(Lq,ΔLf,indices_col)
            
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
                        indices_row,
                        indices_col,
                        ΔLf,
                        syndrome,
                        nothing,
                        nothing,
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
                        indices_row,
                        indices_col,
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
                        indices_row,
                        indices_col,
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
                        indices_row,
                        indices_col,
                        ΔLf,
                        syndrome,
                        sn,
                        nothing,
                        nothing
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

    return log10.(FER), log10.(BER), iters

end