function log_MAP_estimator(N,Lr,Lf, indices_m)

    d = zeros(Int,N)

    for n=1:N
        Ld0 = Lf[n,1]
        Ld1 = Lf[n,2]
        for m in indices_m[n]
            Ld0 += Lr[m,n,1]
            Ld1 += Lr[m,n,2]
        end
        if Ld0 > Ld1
            d[n] = 1
        end
    end
    return d
end