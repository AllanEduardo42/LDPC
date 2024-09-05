################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to estimate the LPCD performance (FER x SNR)

function performance_estimation(c::Vector{Int64},
                                u::Vector{Int64},
                                σ::Vector{Float64},
                                M::Int64, 
                                N::Int64, 
                                indices_n::Vector{Vector{Int64}},
                                indices_m::Vector{Vector{Int64}}, 
                                phi::Vector{Float64},
                                mode::String;
                                nreals = NREALS)

    fer = zeros(length(σ))

    ber = zeros(length(σ),MAX)

    iters = zeros(Int, length(σ), NREALS)

    ΔLf = Vector{Float64}(undef,N)

    ΔLf2 = Vector{Float64}(undef,N)

    Lq = zeros(N,M)

    Lr = zeros(M,N)

    t = Vector{Float64}(undef, N)

    d = Vector{Int64}(undef, N)

    syndrome = Vector{Int64}(undef, M)

    noise = Vector{Float64}(undef, N)

    Lrn = zeros(N)

    sn = ones(Int,N)

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
                @inbounds ΔLf2[i] = -2*SIZE_per_RANGE*t[i]/s
            end

            i = 0
            S = -1
            if mode == "TNH"
                # tanh SPA
                llr_init_q(Lq, N, ΔLf, indices_m)
                S, i = SPA(S,i,M,N,Lr,Lq,indices_n,indices_m,d,ΔLf,syndrome)
            elseif mode == "APP"
                # approximate SPA
                llr_init_q(Lq, N, ΔLf, indices_m)
                S, i = SPA(S,i,M,N,Lr,Lq,indices_n,indices_m,d,ΔLf,syndrome,sn)
            elseif mode == "ALT"
                # alternative SPA
                llr_init_q(Lq, N, ΔLf, indices_m)
                S, i = SPA(S,i,M,N,Lr,Lq,indices_n,indices_m,d,ΔLf,syndrome,sn,Lrn)
            elseif mode == "TAB"
                # lookup-table SPA
                llr_init_q(Lq, N, ΔLf2, indices_m)
                S, i = SPA(S,i,M,N,Lr,Lq,indices_n,indices_m,d,ΔLf2,syndrome,sn,Lrn,phi)
            else
                println("ERROR")
            end
            @inbounds iters[k,j] = i
            if S != 0
                @inbounds fer[k] += 1
            elseif d != c
                @inbounds fer[k] += 1
            end
        end

        @inbounds fer[k] /= NREALS
    end

    return log10.(fer), iters

end