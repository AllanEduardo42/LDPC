# function llr_delta_q(Lq,M,N,indices_m)

#     Lδq = zeros(N,M)
#     S = zeros(N,M)

#     for n=1:N
#         for m in indices_m[n]
#             if Lq[n,m,1] <= Lq[n,m,2]
#                 min_Lq = Lq[n,m,1]
#                 max_Lq = Lq[n,m,2]
#                 S[n,m] = 0
#             else
#                 min_Lq = Lq[n,m,2]
#                 max_Lq = Lq[n,m,1]
#                 S[n,m] = 1
#             end
#             Lδq[n,m] = min_Lq + abs(log(1 - exp(-(max_Lq-min_Lq))))
#             # Lδq[n,m] = min_Lq + get_f_minus(max_Lq-min_Lq)
#             # Lδq[n,m] = min_Lq + f_minus[max_Lq-min_Lq+1]
#         end
#     end

#     return Lδq, S

# end

function llr_simple_horizontal_update(M, N, Lq, indices_n, indices_m)

    Lr = zeros(M,N)

    for m=1:M
        for n in indices_n[m]
            pLr = 1
            for nn in indices_n[m]
                if nn != n
                    pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            Lr[m,n] = 2*atanh(pLr)
        end
    end

    return Lr

end