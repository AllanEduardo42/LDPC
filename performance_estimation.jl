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


    M = length(indices_row)
    N = length(indices_col)

    ############################# preallocation ################################

    fer = zeros(length(σ))

    BER = zeros(MAX,length(σ))

    ber = zeros(MAX)

    iters = zeros(Int, length(σ), NREALS)

    ΔLf = Vector{Float64}(undef,N)

    Lq = H'*0.0

    Lr = H*0.0

    t = Vector{Float64}(undef,N)

    d = Vector{Bool}(undef,N)

    syndrome = Vector{Bool}(undef,M)

    noise = Vector{Float64}(undef,N)

    Lrn = zeros(N)

    sn = ones(Bool,N)

    bit_error = Vector{Bool}(undef,N)

    divisor = NREALS * N

    # BPKS 

    u = 2*c .- 1

    for k in eachindex(σ)

        Random.seed!(SEED)

        s = σ[k]^2
    
        for j=1:nreals

            randn!(noise)

            noise .*= σ[k]

            t .= u .+ noise

            for i in eachindex(t)
                if mode == "TAB"
                    @inbounds ΔLf[i] = -2*SIZE_per_RANGE*t[i]/s
                else
                    @inbounds ΔLf[i] = -2t[i]/s
                end
            end

            llr_init_q!(Lq,ΔLf,indices_col)
            
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
            elseif mode == "APP"
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
            else
                throw(
                    ArgumentError(
                        "$mode is not a valid mode"
                    )
                )
            end
            @fastmath ber ./= divisor
            @inbounds  @fastmath @. BER[:,k] += ber
            @inbounds iters[k,j] = i
            if ~DECODED
                @inbounds fer[k] += 1
            end
        end

        @inbounds @fastmath fer[k] /= NREALS
    end

    return log10.(fer), log10.(BER), iters, Lr

end