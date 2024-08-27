function vertical_update_and_MAP(N,r,q,f,indices_m)

    d = zeros(Int, N)

    for n = 1:N
        d0 = f[n,1]
        d1 = f[n,2]
        for m in indices_m[n]
            d0 *= r[m,n,1]
            d1 *= r[m,n,2]
        end
        if d1 > d0
            d[n] = 1
        end
        for m in indices_m[n]
            q0 = d0 / r[m,n,1]
            q1 = d1 / r[m,n,2]
            α = q0 + q1
            q[n,m,1] = q0/α
            q[n,m,2] = q1/α
        end

    end  

    return q, d

end

function init_q(q,N,f,indices_m)

    for n=1:N
        for m in indices_m[n]
            q[n,m,1] = f[n,1]
            q[n,m,2] = f[n,2]
        end
    end

    return q

end