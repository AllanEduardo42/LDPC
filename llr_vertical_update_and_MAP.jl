function llr_vertical_update_and_MAP(N,Lr,Lq,ΔLf,indices_m)

    d = zeros(Int,N)

    for n = 1:N
        LLd = ΔLf[n]
        for m in indices_m[n]
            LLd += Lr[m,n]
        end
        if LLd < 0
            d[n] = 1
        end
        Lq[n,:] .= LLd
        for m in indices_m[n]
            Lq[n,m] -= Lr[m,n]
        end
    end  

    return Lq, d

end

function llr_init_q(M,N,ΔLf,indices_m)

    Lq = zeros(N,M)

    for n=1:N
        for m in indices_m[n]
            Lq[n,m] = ΔLf[n]
        end
    end

    return Lq

end