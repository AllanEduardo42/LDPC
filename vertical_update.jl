function vertical_update(N,r,q,f,indices_m)

    for n = 1:N
        for m in indices_m[n]
            q0 = f[n,1]
            q1 = f[n,2]
            for mm in indices_m[n]
                if mm != m
                    q0 *= r[mm,n,1]
                    q1 *= r[mm,n,2]
                end
            end
            α = q0 + q1
            q[n,m,1] = q0/α
            q[n,m,2] = q1/α
        end
    end  

    return q

end

function init_q(M,N,f,indices_m)

    q = zeros(N,M,2)

    for n=1:N
        for m in indices_m[n]
            q[n,m,1] = f[n,1]
            q[n,m,2] = f[n,2]
        end
    end

    return q

end