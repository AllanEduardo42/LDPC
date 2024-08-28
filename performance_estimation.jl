################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER x SNR)

function performance_estimation(c::Vector{Int64}, u::Vector{Int64}, σ::Vector{Float64},
                                M::Int64, N::Int64, indices_n::Vector{Vector{Int64}},
                                indices_m::Vector{Vector{Int64}}, h::Matrix{Int64};
                                nreals = NREALS)
    
    Random.seed!(SEED)

    fer = zeros(length(σ))
    iters = zeros(Int, length(σ), NREALS)

    ΔLf = Vector{Float64}(undef,N)

    Lq = zeros(N,M)

    Lr = zeros(M,N)

    t = Vector{Float64}(undef, N)

    d = Vector{Int64}(undef, N)

    pre_syndrome = Vector{Int64}(undef, M)

    syndrome = Vector{Int64}(undef, M)

    noise = Vector{Float64}(undef, N)

    Lrn = zeros(N)

    sn = ones(Int,N)

    for k in eachindex(σ)

        s = σ[k]^2
    
        for j=1:nreals

            randn!(noise)

            noise .*= σ[k]

            t .= u .+ noise

            for i in eachindex(t)
                ΔLf[i] = -2t[i]/s
            end

            llr_init_q(Lq, N, ΔLf, indices_m)

            i = 0
            S = -1

            while S != 0 && i < MAX

                i += 1

                Lr = llr_horizontal_update_alt(Lr,M,Lq,indices_n, Lrn, sn)
                # Lr = llr_horizontal_update(Lr,M,Lq,indices_n)
                Lq, d = llr_vertical_update_and_MAP(Lq, Lr, d, N, ΔLf, indices_m)

                mul!(pre_syndrome, h, d)

                syndrome .= pre_syndrome .% 2 
                
                S = sum(syndrome)

            end

            iters[k,j] = i

            if S != 0
                fer[k] += 1
            elseif d != c
                fer[k] += 1
            end
        end

        fer[k] /= NREALS
    end

    return fer, iters, Lq, Lr

end

function performance_estimation_table(c::Vector{Int64}, u::Vector{Int64}, σ::Vector{Float64},
                                    M::Int64, N::Int64, indices_n::Vector{Vector{Int64}},
                                    indices_m::Vector{Vector{Int64}}, h::Matrix{Int64},
                                    phi::Vector{Float64}; nreals = NREALS)

    Random.seed!(SEED)

    fer = zeros(length(σ))
    iters = zeros(Int, length(σ), NREALS)

    ΔLf = Vector{Float64}(undef,N)

    Lq = zeros(N,M)

    Lr = zeros(M,N)

    t = Vector{Float64}(undef, N)

    d = Vector{Int64}(undef, N)

    pre_syndrome = Vector{Int64}(undef, M)

    syndrome = Vector{Int64}(undef, M)

    noise = Vector{Float64}(undef, N)

    Lrn = zeros(N)

    sn = ones(Int,N)

    for k in eachindex(σ)

        s = σ[k]^2

        for j=1:nreals

            randn!(noise)

            noise .*= σ[k]

            t .= u .+ noise

            for i in eachindex(t)
                ΔLf[i] = -2*SIZE_per_RANGE*t[i]/s
            end

            llr_init_q(Lq, N, ΔLf, indices_m)

            i = 0
            S = -1

            while S != 0 && i < MAX

                i += 1

                Lr =  llr_horizontal_update_table(Lr,M,Lq,indices_n, Lrn, sn, phi)
                Lq, d = llr_vertical_update_and_MAP(Lq, Lr, d, N, ΔLf, indices_m)

                pre_syndrome .*= 0
                for m=1:M
                    for n in indices_n[m]
                        pre_syndrome[m] += d[n]
                    end
                end

                syndrome .= pre_syndrome .% 2 

                S = sum(syndrome)

            end

            iters[k,j] = i

            if S != 0
                fer[k] += 1
            elseif d != c
                fer[k] += 1
            end
        end

        fer[k] /= NREALS
    end

    return fer, iters, Lq, Lr

end