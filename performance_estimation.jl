################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER x SNR)

function performance_estimation(c::Vector{Int64},
                                u::Vector{Int64},
                                σ::Vector{Float64},
                                H::Matrix{Int64},
                                indices_n::Vector{Vector{Int64}},
                                indices_m::Vector{Vector{Int64}}, 
                                phi::Vector{Float64},
                                mode::String;
                                nreals = NREALS)


    M = length(indices_n)
    N = length(indices_m)

    fer = zeros(length(σ))

    BER = zeros(MAX,length(σ))

    ber = zeros(MAX)

    iters = zeros(Int, length(σ), NREALS)

    ΔLf = Vector{Float64}(undef,N)

    Lq = sparse(H'*1e-16)

    Lr = sparse(H*1e-16)

    # Lq = H'*0.0

    # Lr = H*0.0

    t = Vector{Float64}(undef, N)

    d = Vector{Int64}(undef, N)

    syndrome = Vector{Int64}(undef, M)

    noise = Vector{Float64}(undef, N)

    Lrn = zeros(N)

    sn = ones(Int,N)

    x = BitArray(undef, length(c))

    divisor = NREALS * N

    for k in eachindex(σ)

        Random.seed!(SEED)

        s = σ[k]^2
    
        for j=1:nreals

            randn!(noise)

            noise .*= σ[k]

            t .= u .+ noise

            for i in eachindex(t)
                @inbounds ΔLf[i] = -2t[i]/s
            end

            for i in eachindex(t)
                if mode == "TAB"
                    @inbounds ΔLf[i] = -2*SIZE_per_RANGE*t[i]/s
                else
                    @inbounds ΔLf[i] = -2t[i]/s
                end
            end

            llr_init_q!(Lq,ΔLf,indices_m)
            
            i = 0
            DECODED = false
            if mode == "TNH"
                # tanh SPA
                DECODED, i = SPA!(d,ber,c,x,Lr,Lq,indices_n,indices_m,ΔLf,
                                  syndrome,nothing,nothing,nothing)
            elseif mode == "APP"
                # approximate SPA
                DECODED, i = SPA!(d,ber,c,x,Lr,Lq,indices_n,indices_m,ΔLf,
                                 syndrome,sn,nothing,nothing)
            elseif mode == "ALT"
                # alternative SPA
                DECODED, i = SPA!(d,ber,c,x,Lr,Lq,indices_n,indices_m,ΔLf,
                                  syndrome,sn,Lrn,nothing)
            elseif mode == "TAB"
                # lookup-table SPA
                DECODED, i = SPA!(d,ber,c,x,Lr,Lq,indices_n,indices_m,ΔLf,
                                  syndrome,sn,Lrn,phi)
            else
                throw(ArgumentError("invalid mode"))
            end
            ber ./= divisor
            @inbounds @. BER[:,k] += ber
            @inbounds iters[k,j] = i
            if ~DECODED
                @inbounds fer[k] += 1
            end
        end

        @inbounds fer[k] /= NREALS
    end

    return log10.(fer), log10.(BER), iters

end