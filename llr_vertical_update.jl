function llr_vertical_update(N,Lr,Lq,Lf,indices_m)

    for n = 1:N
        for m in indices_m[n]
            Lq[n,m] = Lf[n,1] - Lf[n,2]
            for mm in indices_m[n]
                if mm != m
                    Lq[n,m] += Lr[mm,n]
                end
            end
        end
    end  

    return Lq

end

function llr_init_q(M,N,Lf,indices_m)

    Lq = zeros(N,M)

    for n=1:N
        for m in indices_m[n]
            Lq[n,m] = Lf[n,1] - Lf[n,2]
        end
    end

    return Lq

end