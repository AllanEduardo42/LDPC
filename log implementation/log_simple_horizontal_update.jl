function log_delta_q(Lq,M,N,indices_m)

    Lδq = zeros(N,M)
    S = zeros(N,M)

    for n=1:N
        for m in indices_m[n]
            if Lq[n,m,1] <= Lq[n,m,2]
                min_Lq = Lq[n,m,1]
                max_Lq = Lq[n,m,2]
                S[n,m] = 0
            else
                min_Lq = Lq[n,m,2]
                max_Lq = Lq[n,m,1]
                S[n,m] = 1
            end
            Lδq[n,m] = min_Lq + abs(log(1 - exp(-(max_Lq-min_Lq))))
            # Lδq[n,m] = min_Lq + get_f_minus(max_Lq-min_Lq)
            # Lδq[n,m] = min_Lq + f_minus[max_Lq-min_Lq+1]
        end
    end

    return Lδq, S

end

function log_simple_horizontal_update(M, N, Lq, indices_n, indices_m)

    Lr = zeros(M,N,2)

    Lδq, S = log_delta_q(Lq,M,N,indices_m)

    for m=1:M
        for n in indices_n[m]
            Lδr = 0
            sδr = 0
            for nn in indices_n[m]
                if nn != n
                    Lδr += Lδq[nn,m]
                    sδr += S[nn,m]
                end
            end
            if iseven(sδr)
                Lr[m,n,1] = log(2) - abs(log(1 + exp(-Lδr)))
                Lr[m,n,2] = log(2) + abs(log(1 - exp(-Lδr)))

                # Lr[m,n,1] = SCOPE_LN2 - get_f_plus(Lδr)
                # Lr[m,n,2] = SCOPE_LN2 + get_f_minus(Lδr)

                # Lr[m,n,1] = SCOPE_LN2 - f_plus[Lδr+1]
                # Lr[m,n,2] = SCOPE_LN2 + f_minus[Lδr+1]
            else
                Lr[m,n,1] = log(2) + abs(log(1 - exp(-Lδr)))
                Lr[m,n,2] = log(2) - abs(log(1 + exp(-Lδr)))

                # Lr[m,n,1] = SCOPE_LN2 + get_f_minus(Lδr)
                # Lr[m,n,2] = SCOPE_LN2 - get_f_plus(Lδr)

                # Lr[m,n,1] = SCOPE_LN2 + f_minus[Lδr+1]
                # Lr[m,n,2] = SCOPE_LN2 - f_plus[Lδr+1]
            end
        end
    end

    return Lr

end