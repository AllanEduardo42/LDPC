function log_vertical_update(N,Lr,Lq,Lf,indices_m)

    for n = 1:N
        for m in indices_m[n]
            Lq0 = Lf[n,1]
            Lq1 = Lf[n,2]
            for mm in indices_m[n]
                if mm != m
                    Lq0 += Lr[mm,n,1]
                    Lq1 += Lr[mm,n,2]
                end
            end
            if Lq0 >= Lq1
                min_Lq = Lq1
                max_Lq = Lq0
            else
                min_Lq = Lq0
                max_Lq = Lq1
            end
            α = min_Lq - abs(log(1 + abs(exp(-(max_Lq-min_Lq)))))
            # α = min_Lq - get_f_plus(max_Lq-min_Lq)
            # α = min_Lq - f_plus[max_Lq-min_Lq+1]
            Lq[n,m,1] = Lq0 - α
            Lq[n,m,2] = Lq1 - α
        end
    end  

    return Lq

end

function log_init_q(M,N,Lf,indices_m)

    Lq = zeros(N,M,2)

    for n=1:N
        for m in indices_m[n]
            Lq[n,m,1] = Lf[n,1]
            Lq[n,m,2] = Lf[n,2]
        end
    end

    return Lq

end