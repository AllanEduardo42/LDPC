function findindices_M(H,N)

    indices_m = Vector{Vector}(undef, N)
    for n=1:N
        indices_m[n] = findall(x -> x == 1, H[:,n])
    end

    return indices_m

end

function findindices_N(H,M)

    indices_n = Vector{Vector}(undef, M)
    for m=1:M
        indices_n[m] = findall(x -> x == 1, H[m,:])
    end

    return indices_n

end

function normalize(f, N)

    for n=1:N
        α = f[n,1] + f[n,2]
        f[n,1] = f[n,1]/α
        f[n,2] = f[n,2]/α
    end

end