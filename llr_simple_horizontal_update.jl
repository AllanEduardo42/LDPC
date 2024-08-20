function llr_simple_horizontal_update(M, N, Lq, indices_n)

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