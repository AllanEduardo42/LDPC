function llr_MAP_estimator(N,Lr,Lf, indices_m)

    d = zeros(Int,N)

    for n=1:N
        Ld = Lf[n,1] - Lf[n,2]
        for m in indices_m[n]
            Ld += Lr[m,n]
        end
        if Ld < 0
            d[n] = 1
        end
    end
    return d
end