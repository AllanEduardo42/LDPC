function performance_estimation(r, σ, M, N, indices_n, indices_m, h; nreals = NREALS)

    fer = zeros(length(σ))
    iters = zeros(Int, length(σ), NREALS)

    # f = Matrix{Float64}(undef,N,2)

    # Lf = Matrix{Float64}(undef,N,2)

    ΔLf = Vector{Float64}(undef,N)

    Lq = zeros(N,M)

    Lr = zeros(M,N)

    t = Vector{Float64}(undef, N)

    d = Vector{Int64}(undef, N)

    pre_syndrome = Vector{Int64}(undef, M)

    syndrome = Vector{Int64}(undef, M)

    noise = Vector{Float64}(undef, N)

    for k in eachindex(σ)

        c = 1/(sqrt(2π)*σ[k])
        s = 2σ[k]^2
    
        for j=1:nreals

            # Lq .= zeros(N,M)
            # Lr .= zeros(M,N)
            # d .= zeros(Int,N)

            randn!(noise)

            noise .*= σ[k]

            # t .= r .+ noise

            t = [1.3129, 2.6584, 0.7413, 2.1745, 0.5981, −0.8323, −0.3962, −1.7586,
            1.4905, 0.4084, −0.9290, 1.0765]

            for i in eachindex(t)
                # f[i,1] = c*exp(-(t[i]+1)^2/s)
                # f[i,2] = c*exp(-(t[i]-1)^2/s)
                ΔLf[i] = -4t[i]/s
            end

            # normalize(f,N)
            # Lf .= log.(f)
            # ΔLf = view(Lf,:,1) - view(Lf,:,2)

            llr_init_q(Lq, N, ΔLf, indices_m)

            i = 0
            S = -1

            while S != 0 && i < MAX

                i += 1

                Lr = llr_simple_horizontal_update(Lr,M,Lq,indices_n)
                Lq, d = llr_vertical_update_and_MAP(Lq, Lr, d, N, ΔLf, indices_m)

                println("d=",d)

                mul!(pre_syndrome, h, d)

                syndrome .= rem.(pre_syndrome,2)    
                
                S = sum(syndrome)

            end

            iters[k,j] = i

            if S !== 0
                fer[k] += 1
            end
        end

        fer[k] /= NREALS
    end

    return fer, iters, Lq, Lr

end