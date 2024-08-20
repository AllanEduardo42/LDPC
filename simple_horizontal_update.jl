function simple_horizontal_update(M, N, q, indices_n)

    r = zeros(M,N,2)

    δq = q[:,:,1]-q[:,:,2]

    for m=1:M
        for n in indices_n[m]
            δr = 1
            for nn in indices_n[m]
                if nn != n
                    δr *= δq[nn,m]
                end
            end
            r[m,n,1] = 0.5*(1+δr)
            r[m,n,2] = 0.5*(1-δr)
        end
    end

    return r

end