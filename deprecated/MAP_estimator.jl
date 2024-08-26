function MAP_estimator(N,r,f, indices_m)

    d = zeros(Int,N)

    for n=1:N
        d0 = f[n,1]
        d1 = f[n,2]
        for m in indices_m[n]
            d0 *= r[m,n,1]
            d1 *= r[m,n,2]
        end
        if d1 > d0
            d[n] = 1
        end
    end
    return d
end